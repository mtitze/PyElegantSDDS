from .main import ElegantRunToolkit

class radiation(ElegantRunToolkit):

    def add_radiation_damping(self, isr=1, synch_rad=1, isr1part=1, use_rad_dist=1, use_bend_opening_angle=1, **kwargs):
        """Add radiation damping for CSBEND, KQUAD and KSEXT.
        
        isr: also set ISR (include incoherent synchrotron radiation) to the given value.
        """
        self.er.add_alter_elements(name="*", type="CSBEND", item="SYNCH_RAD", value=synch_rad, allow_missing_elements=1)
        self.er.add_alter_elements(name="*", type="CSBEND", item="USE_RAD_DIST", value=use_rad_dist, allow_missing_elements=1)
        if use_rad_dist == 1:
            if use_bend_opening_angle == 1:
                self.er.add_alter_elements(name="*", type="CSBEND", item="ADD_OPENING_ANGLE", value=use_bend_opening_angle, allow_missing_elements=1)
            else:
                print ('use_bend_opening_angle = 1 only possible if use_rad_dist = 1') # see https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu126.html#x137-13600010.19)
        self.er.add_alter_elements(name="*", type="CSBEND", item="ISR", value=isr, allow_missing_elements=1)
        
        self.er.add_alter_elements(name="*", type="KQUAD", item="SYNCH_RAD", value=synch_rad, allow_missing_elements=1)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR", value=isr, allow_missing_elements=1)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR1PART", value=isr1part, allow_missing_elements=1)
        
        self.er.add_alter_elements(name="*", type="KSEXT", item="SYNCH_RAD", value=synch_rad, allow_missing_elements=1)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR", value=isr, allow_missing_elements=1)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR1PART", value=isr1part, allow_missing_elements=1)