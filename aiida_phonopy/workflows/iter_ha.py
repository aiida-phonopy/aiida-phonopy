import sys

class IterHarmonicApproximation(WorkChain):
    _next_workchain_string = 'phonopy.phonon'
    _next_workchain = WorkflowFactory(_next_workchain_string)

    @classmethod
    def define(cls, spec):
        super(IterHarmonicApproximation, cls).define(spec)
        spec.expose_inputs(cls._next_workchain)
        spec.outline(
            cls.initialize,
        )

    def initialize(self):
        self.report("initialize")
