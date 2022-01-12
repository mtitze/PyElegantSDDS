from .main import ElegantRunToolkit

class radiation(ElegantRunToolkit):

    def add_radiation_damping(self, isr=1, synch_rad=1, isr1part=1, use_rad_dist=1, **kwargs):
        """Add radiation damping for CSBEND, KQUAD and KSEXT.
        
        isr: also set ISR (include incoherent synchrotron radiation) to the given value.
        """
        self.er.add_alter_elements(name="*", type="CSBEND", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="CSBEND", item="USE_RAD_DIST", value=use_rad_dist)
        self.er.add_alter_elements(name="*", type="CSBEND", item="ISR", value=isr)
        
        self.er.add_alter_elements(name="*", type="KQUAD", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR", value=isr)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR1PART", value=isr1part)
        
        self.er.add_alter_elements(name="*", type="KSEXT", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR", value=isr)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR1PART", value=isr1part)