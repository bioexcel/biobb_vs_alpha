#!/usr/bin/env python

"""Python wrapper module for the GROMACS genrestr module
"""
import sys
import re
import json
import os
import configuration.settings as settings
from command_wrapper import cmd_wrapper
from tools import file_utils as fu


class Genrestr(object):
    """Wrapper class for the GROMACS genrestr module.

    Args:
        input_structure_path (str): Path to the input structure PDB/GRO/TPR file.
        input_ndx_path (str): Path to the input index NDX file.
        input_top_zip_path (str): Path the input TOP topology in zip format.
        output_top_zip_path (str): Path the output TOP topology in zip format.
        properties (dic):
            | **output_top_path** (*str*): Path the output TOP file.
            | **output_itp_path** (*str*): Path to the output include for topology ITP file.
            | **force_constants** (*float[3]*): Array of three floats defining the force constants
    """

    def __init__(self, input_structure_path, input_ndx_path, input_top_zip_path,
                 output_top_zip_path, properties, **kwargs):
        if isinstance(properties, basestring):
            properties=json.loads(properties)
        self.input_structure_path = input_structure_path
        self.input_ndx_path = input_ndx_path
        self.input_top_zip_path = input_top_zip_path
        self.output_top_zip_path = output_top_zip_path
        self.output_itp_path = properties.get('output_itp_path','restrain.itp')
        self.output_top_path = properties.get('output_top_path','restrain.top')
        self.force_constants = properties.get('force_constants','500 500 500')
        self.restricted_group = properties.get('restricted_group', 'system')
        self.gmx_path = properties.get('gmx_path',None)
        self.mutation = properties.get('mutation',None)
        self.step = properties.get('step',None)
        self.path = properties.get('path','')
        self.mpirun = properties.get('mpirun',False)
        self.mpirun_np = properties.get('mpirun_np',None)

    def launch(self):
        """Launches the execution of the GROMACS pdb2gmx module.
        """
        out_log, err_log = fu.get_logs(path=self.path, mutation=self.mutation, step=self.step)
        self.output_top_path = fu.add_step_mutation_path_to_name(self.output_top_path, self.step, self.mutation)
        self.output_itp_path = fu.add_step_mutation_path_to_name(self.output_itp_path, self.step, self.mutation)

        gmx = "gmx" if self.gmx_path is None else self.gmx_path
        cmd = [gmx, "genrestr", "-f", self.input_structure_path,
               "-n", self.input_ndx_path, "-o", self.output_itp_path,
               "-fc", self.force_constants]

        if self.mpirun_np is not None:
            cmd.insert(0, str(self.mpirun_np))
            cmd.insert(0, '-np')
        if self.mpirun:
            cmd.insert(0, 'mpirun')

        if self.mpirun:
            cmd.append('<<<')
            cmd.append('\"'+self.restricted_group+'\"')
        else:
            cmd.insert(0, '|')
            cmd.insert(0, '\"'+self.restricted_group+'\"')
            cmd.insert(0, 'echo')

        command = cmd_wrapper.CmdWrapper(cmd, out_log, err_log)
        returncode = command.launch()

        fu.unzip_top(zip_file=self.input_top_zip_path, top_file=self.output_top_path)
        out_log.info('Unzip: '+ self.input_top_zip_path + ' to: '+self.output_top_path)
        with open(self.output_top_path, 'r') as fin:
            for line in fin:
                if line.startswith('#ifdef POSRES'):
                    itp_name = re.findall('"([^"]*)"',fin.next())[0]
                    out_log.debug('itp_name: '+itp_name)
                    break

        # with open(self.output_top_path, 'r') as fin:
        #     data = fin.read().splitlines(True)
        #     index = data.index('#ifdef POSRES\n')
        #     data[index+2] = 'system\n'
        #     data.insert(index, '\n')
        #     data.insert(index, '#endif\n')
        #     data.insert(index, '#include "'+self.output_itp_path+'"\n')
        #     data.insert(index, '#ifdef CUSTOM_POSRES\n')
        #     data.insert(index, '; Include Position restraint file\n')
        #     # data.insert(index, '#include "'+self.output_itp_path+'"\n')
        #     # data.insert(index, '; Include genrestr generated itp\n')
        # with open(self.output_top_path, 'w') as fout:
        #     fout.writelines(data)

        with open(self.output_itp_path, 'r') as fin:
            data = fin.read().splitlines(True)
            # data.insert(0, '\n')
            # data.insert(0, 'system    3\n')
            # data.insert(0, ';Name    nrexcl\n')
            # data.insert(0, '[ system ]\n')
        with open(itp_name, 'w') as fout:
            fout.writelines(data)
        os.remove(self.output_itp_path)


        # zip topology
        fu.zip_top(self.output_top_path, self.output_top_zip_path, remove_files=False)
        out_log.info('Zip: '+ self.output_top_path +' to: '+ self.output_top_zip_path)

        return returncode
#Creating a main function to be compatible with CWL
def main():
    system=sys.argv[1]
    step=sys.argv[2]
    properties_file=sys.argv[3]
    prop = settings.YamlReader(properties_file, system).get_prop_dic()[step]
    Pdb2gmx(input_structure_pdb_path=sys.argv[4],
            output_gro_path=sys.argv[5],
            output_top_zip_path=sys.argv[6],
            properties=prop).launch()

if __name__ == '__main__':
    main()
