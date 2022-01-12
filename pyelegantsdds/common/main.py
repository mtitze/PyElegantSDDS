from pyelegantsdds.elegantrun import ElegantRun
import copy

class ElegantRunToolkit:
    '''
    Class to perform advanced operations using the ElegantRun functionality.
    '''
    def __init__(self, er: ElegantRun):
        self.er = er

    def copy_elegant_run(self, new_rootname=''):
        # create new instance of current elegant run, using the same init parameters
        # as the original one, with a new rootname.
        if len(new_rootname) == 0:
            new_rootname = self.er.rootname + '_copy'
        return ElegantRun(sif=self.er.sif, lattice=self.er.lattice, energy_gev=self.er.energy_gev,
                             beamline=self.er.beamline, rootname=new_rootname,
                             parallel=self.er.parallel)
        
    def copy(self, new_rootname=''):
        # create a copy of the current elegant run, including its command list
        er_cp = self.copy_elegant_run(new_rootname)
        er_cp.commandfile.commandlist = copy.deepcopy(self.er.commandfile.commandlist)
        er_cp.commandfile.history = copy.deepcopy(self.er.commandfile.history)
        return er_cp
