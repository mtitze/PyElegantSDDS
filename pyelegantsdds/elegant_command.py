# -*- coding: utf-8 -*-

"""
Module pyelegantsdds.elegant_command 
=================================================================

A module containing the ElegantCommandFile class.

"""

from itertools import count


class ElegantCommandFile:
    _COMMANDLIST = (
        "alter_elements",
        "amplification_factors",
        "analyze_map",
        "aperture_data",
        "bunched_beam",
        "change_particle",
        "chromaticity",
        "closed_orbit",
        "correct",
        "correction_matrix_output",
        "correct_tunes",
        "coupled_twiss_output",
        "divide_elements",
        "error_element",
        "error_control",
        "find_aperture",
        "floor_coordinates",
        "frequency_map",
        "global_settings",
        "insert_elements",
        "insert_sceffects",
        "linear_chromatic_tracking_setup",
        "link_control",
        "link_elements",
        "load_parameters",
        "matrix_output",
        "modulate_elements",
        "moments_output",
        "momentum_aperture",
        "optimize",
        "optimization_constraint",
        "optimization_covariable",
        "optimization_setup",
        "parallel_optimization_setup",
        "optimization_term",
        "optimization_variable",
        "print_dictionary",
        "ramp_elements",
        "rf_setup",
        "replace_elements",
        "rpn_expression",
        "rpn_load",
        "run_control",
        "run_setup",
        "sasefel",
        "save_lattice",
        "sdds_beam",
        "semaphores",
        "slice_analysis",
        "subprocess",
        "steering_element",
        "touschek_scatter",
        "transmute_elements",
        "twiss_analysis",
        "twiss_output",
        "track",
        "tune_shift_with_amplitude",
        "vary_element",
    )

    def __init__(self, filename):
        self.commandlist = []
        self.filename = filename
        self.history = []

    def checkCommand(self, typename):
        """
        Check if a command is a valid
        Elegant command.

        Parameters:
        ----------
        typename    : str
                command type

        """
        return typename.lower() in self._COMMANDLIST
    
    def exists(self, command):
        '''
        Check if a command is present in the current command list.
        
        Parameters
        ----------
        command: str
            Name of the command.
        '''
        cmd = command.lower()
        if not cmd in [cmd['NAME'] for cmd in self.commandlist]:
            raise RuntimeError(f'Command {cmd} required.')

    def addCommand(self, command, note="", **params):
        """
        Method to add a command to the command file.
        Generates a dict to reconstruct command string,
        that can be added to the command file.
        
        values of the form 'None' or '""' are ignored.

        Parameters:
        ----------
        command     : str
                valid Elegant command
        note        : str
                extra info
        """
        # check if it is a valid Elegant command
        if not self.checkCommand(command):
            raise RuntimeError("Command '{}' not recognized.".format(command))

        ignored_values = ["", None]
        
        # init command dict
        thiscom = {}

        thiscom["NAME"] = command.lower()
        thiscom["NOTE"] = note

        # add command parameters to the command dict
        for k, v in params.items():
            if v in ignored_values:
                continue
            thiscom[k] = v

        # add the command dict to the command list
        self.commandlist.append(thiscom)

    def modifyCommand(self, commandname, mode="last", **params):
        """
        Method to modify a command already in the command list.
        If the command is used multiple times, the command you
        want to change can be select by using the strings "first" or "last"
        or by giving the integer index of the command.

        Parameters:
        ----------
        commandname     : str
                name of the command to modify
        mode            : str | int, default last
                if multiple commands select which one to change
        params          : dict
                command parameters - updated

        """
        # create command index list to be able
        # to choose which one to change
        indlist = []
        i = 0

        for command in self.commandlist:
            if command["NAME"] == commandname.lower():
                indlist.append(i)
            i += 1

        # find command to updated given in mode
        # allowed are : first, last or integer
        if len(indlist) == 1:
            for k, v in params.items():
                self.commandlist[indlist[0]][k] = v
            return
        elif len(indlist) == 0:
            raise RuntimeError("Command '{}' not found".format(commandname))
        else:
            if mode == "last":
                for k, v in params.items():
                    self.commandlist[indlist[-1]][k] = v
                return
            elif mode == "first":
                for k, v in params.items():
                    self.commandlist[indlist[0]][k] = v
                return
            elif isinstance(mode, int):
                for k, v in params.items():
                    self.commandlist[indlist[mode]][k] = v
                return
            else:
                raise RuntimeError("The mode is invalid")

    def repeatCommand(self, commandname, mode="last"):
        """
        Method to repeat a commend already in the
        command list.


        Parameters:
        ----------
        commandname     : str
                name of the command to modify
        mode            : str | int , default last
                if multiple commands select which one to change

        """
        # create index list for the command to be repeated
        indlist = []
        i = 0
        for command in self.commandlist:
            if command["NAME"] == commandname.lower():
                indlist.append(i)
            i += 1

        if len(indlist) == 1:
            self.commandlist.append(self.commandlist[indlist[0]])
            return
        elif len(indlist) == 0:
            raise RuntimeError("Command '{}' not found".format(commandname))
        else:
            if mode == "last":
                self.commandlist.append(self.commandlist[indlist[-1]])
            elif mode == "first":
                self.commandlist.append(self.commandlist[indlist[0]])
            elif isinstance(mode, int):
                self.commandlist.append(self.commandlist[indlist[mode]])
            else:
                raise RuntimeError("The mode is invalid")

    def clearHistory(self):
        """
        Clear command history.
        """
        self.history = {}

    def addHistory(self, cmdlist=None):
        """
        Add commands to history.

        Parameters:
        -----------
        list | str : cmdlist
                command or list of commands to add to history.
        """
        if cmdlist == None:
            cmdlist = self.commandlist
        if not isinstance(cmdlist, list):
            raise RuntimeError('commandlist not a list.')
        self.history.append(cmdlist)

    def clear(self):
        """Clear command list."""
        # check if commadnlist empty if not add to history
        if len(self.commandlist) > 0:
            self.addHistory(self.commandlist)

        # clear commandlist
        self.commandlist = []

    def remove_command(self, commandname, mode="last"):
        # create command index list to be able
        # to choose which one to change
        clist = [d.get("NAME") for d in self.commandlist]
        indlist = [i for i, x in enumerate(clist) if x == commandname]

        if mode == "last":
            self.commandlist = (
                self.commandlist[: indlist[-1]] + self.commandlist[indlist[-1] + 1 :]
            )
        elif mode == "first":
            self.commandlist = self.commandlist[: indlist[0]] + self.commandlist[indlist[0] + 1 :]
        else:
            self.commandlist = (
                self.commandlist[: indlist[mode]] + self.commandlist[indlist[mode] + 1 :]
            )

    def write(self, outputfilename="", mode="w", verbose=True):
        """
        Method to write the command file to external file.

        Parameters:
        -----------
        outputfilename: str
                filename to write to
        mode : str
                python file write mode
        """
        if outputfilename != "":
            self.filename = outputfilename

        # writing the command file
        with open(self.filename, mode) as outfile:
            # looping over the commands
            for command in range(len(self.commandlist)):
                outfile.write("&{}\n".format(self.commandlist[command]["NAME"]))

                if self.commandlist[command]["NOTE"] != "":
                    outfile.write("! {}\n".format(self.commandlist[command]["NOTE"]))

                for k, v in self.commandlist[command].items():
                    if k != "NAME" and k != "NOTE":
                        if len(k) <= 19:
                            outfile.write("\t{:20s}= {},\n".format(k, v))
                        elif len(k) <= 29:
                            outfile.write("\t{:30s}= {},\n".format(k, v))
                        else:
                            outfile.write("\t{:40s}= {},\n".format(k, v))
                outfile.write("&end\n\n")
        if verbose:
             print (f'Writing command file\n{self.filename}\n')
