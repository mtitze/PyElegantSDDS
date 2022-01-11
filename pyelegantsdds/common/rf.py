from .main import ElegantRunToolkit

from pyelegantsdds.sdds import SDDS
import numpy as np

class rftools(ElegantRunToolkit):
    '''
    Class to perform RF modifications to an Elegant run.
    '''
        
    def add_rf_cavity(self, length, voltage, frequency, phase, **kwargs):
        """Add generic rf cavity (at lattice start by default)."""
        self.er.commandfile.exists('run_setup')
        self.er.commandfile.addCommand(
            "insert_elements",
            name=kwargs.get("name", "RF"),
            type=kwargs.get("type", "RFCA"),
            exclude="",
            s_start=kwargs.get("s_start", -1),
            s_end=kwargs.get("s_end", -1),
            skip=kwargs.get("skip", 1),
            insert_before=kwargs.get("insert_before", 0),
            add_at_end=kwargs.get("add_at_end", 0),
            add_at_start=kwargs.get("add_at_start", 1),
            #element_def=kwargs.get("element_def", r'"RF: RFCA, "')
            element_def=f'"RF: RFCA, L={length}, VOLT={voltage}, FREQ={frequency}, PHASE={phase}"'
            )

    def get_rf_get_freq_and_phase(self, total_voltage: float, harmonic: int, rad: bool=False, **kwargs):
        """Add simulation single turn to get synchronuous freq and phase to be used in next simulation.
        Produces temp.twi and temp.out that can be used with tags.
        Parameters
        ----------
        total_voltage : float, optional
            Total sum of RF voltages, by default 4*375e3
        harmonic : int, optional
            harmonic number, by default 400
        rad : bool, optional
            radiation flag, by default False
        """
        
        # create a copy of the current Elegant run to get the RF frequency and the phase.
        rootname2 = f'{self.er.rootname}_rf_freq_and_phase'
        er2 = self.copy_elegant_run(rootname2)
        er2.add_run_setup(**kwargs)
        if rad:
            rad_er2_tools = radiation(er2)
            rad_er2_tools.add_radiation_damping(**kwargs)

        er2.add_twiss_output()
        er2.add_rf_setup(total_voltage=total_voltage, harmonic=harmonic, track_for_frequency=1)
        er2.commandfile.addCommand("run_control")
        er2.commandfile.addCommand("track")
        er2.run()
        return rootname2
        
    def rf_handling(self, rf_status, total_voltage, one_turn_rootname='', **kwargs):
        
        # rf_status: 0 no RF changes
        # rf_status: 1: insert generic RF cavity
        # rf_status: 2: use/modify existing RF cavitie(s)
        
        if rf_status == 1:
            harmonic = kwargs.get('harmonic')
            rf_length = kwargs.get('rf_length', 0)
            rf_freq = kwargs.get('rf_frequency')
            assert rf_freq != None
            rf_phase = kwargs.get('rf_phase', 180)
            
            if kwargs.get("rad", False) and len(one_turn_rootname) > 0:
                twissfile = f"{one_turn_rootname}.twi"
                print (f'Loading U0 from {twissfile} to modify rf phase ...')
                sdds_one_turn_twi = SDDS(self.er.sif, twissfile, 0, rootname=f'{one_turn_rootname}')
                twiss_params = sdds_one_turn_twi.getParameterValues()
                U0 = twiss_params.U0*1e6 # convert MeV to eV
                phase_correction = -180*np.arcsin(U0/total_voltage)
            else:
                phase_correction = 0
            rf_phase += phase_correction
            print (f'RF phase: {rf_phase} deg')
            
            #self.er.add_twiss_output() # Commented out due to Elegant warning: To avoid possible calculation errors, insert_elements commands should immediately follow run_setup
            self.add_rf_cavity(length=rf_length, voltage=total_voltage, frequency=rf_freq, phase=rf_phase, **kwargs)
            #self.add_rf_setup(filename="%s.rf", **kwargs)
            #self.add_alter_elements(name="*", type="RFCA", item="PHASE", value=f'"({phase})"') # bugfix: need to change PHASE sepeartely, because rf_setup apparently computes the frequency not independently!
        
        if rf_status == 2:
            harmonic = kwargs.get('harmonic')
            n_cavities = kwargs.get('n_cavities')
            assert n_cavities > 0 # may need to get this from the lattice automatically...
            assert harmonic != None

            # get frequency
            if len(one_turn_rootname) == 0:
                print ('one_turn_rootname not provided. Attempting to obtain the frequency from tracking ...')
                one_turn_rootname = self.get_rf_get_freq_and_phase(total_voltage=total_voltage, harmonic=harmonic, rad=False) # run with track_for_frequency=1; TODO: need to discuss case rad=True
            
            self.er.add_twiss_output()
            self.er.commandfile.addCommand("rpn_load", filename=f'{one_turn_rootname}.out', tag="ref", use_row=0) ### filename exists???

            # get energy loss
            self.er.commandfile.addCommand("rpn_load", filename=f"{one_turn_rootname}.twi", tag="reftwi", load_parameters=1)

            # update RF settings
            # VOLTAGE
            self.er.add_alter_elements(name="*", type="RFCA", item="VOLT", value=total_voltage/n_cavities)
            # FREQUENCY
            # N.B. Elegant uses RPN (reverse polnic notation) for formula expressions.
            # The definition of the operator 'rec' (= reciprocal) is hereby given in the file
            # mti/containers/elegant/3rdparty/defns.rpn
            # The following operation sets the frequencies in all RF-cavities according to 1/t_ref*harmonic.
            # !!! Note that in reality, these cavities may operate with different frequencies.
            self.er.add_alter_elements(name="*", type="RFCA", item="FREQ", value=f'"(ref.t rec {harmonic} * )"') # bugfix: brackets and " were in the wrong order!

            # PHASE
            # Note: reftwi.U0 is the energy loss per turn, in Elegant given in MeV
            self.er.add_alter_elements(name="*", type="RFCA", item="PHASE", value=f'"(180 reftwi.U0 1e6 * {total_voltage} / dasin -)"') # N.B. dasin: degree-arcus-sinus: dasin() = 180/pi*arcsin() yields degreees
            # PHASE 2nd varation
            # self.add_alter_elements(name="*", type="RFCA", item="PHASE", value='"(180)"')

            #self.add_twiss_output()