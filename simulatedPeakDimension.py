
class SimulatedPeakDimension(object):
    """Dimension of an expected peak."""
    def __init__(self, simulatedPeak, atom):
        super(SimulatedPeakDimension, self).__init__()
        self.simulatedPeak = simulatedPeak
        self.atom = atom
        self.shift = None
        self.shift_is_assigned = False
        self.update_shift()

    def update_shift(self):
        '''Update shift. Either gets the shift
           from a shift list, if present. Otherwise
           takes average shift in the reference
           database.

        '''
        self.shift = self.get_shift_from_list()
        if self.shift:
            self.shift_is_assigned = True
            return
        self.shift_is_assigned = False
        self.shift = self.get_average_reference_shift()

    def get_shift_from_list(self):
        '''Get the shift for the atom in this peak
           dimension from a shift list.

        '''
        shiftList = self.simulatedPeak.simulated_partial_spectrum.shiftList
        if not shiftList:
            return
        shifts = []
        if self.atom.atomSet and self.atom.atomSet.resonanceSets:
            for resonanceSet in self.atom.atomSet.resonanceSets:
                for resonance in resonanceSet.resonances:
                    shift = resonance.findFirstShift(parentList=shiftList)
                    if shift:
                        shifts.append(shift.value)
        if shifts:
            return sum(shifts)/float(len(shifts))


    def get_average_reference_shift(self):
        '''Get shift from reference data base.'''
        atomName = self.atom.chemAtom.name
        project = self.atom.residue.chain.molSystem.parent
        ccpCode = self.atom.residue.ccpCode
        nmrRefStore = project.findFirstNmrReferenceStore(molType='protein',
                                                         ccpCode=ccpCode)
        chemCompNmrRef = nmrRefStore.findFirstChemCompNmrRef(sourceName='RefDB')
        chemCompVarNmrRef = chemCompNmrRef.findFirstChemCompVarNmrRef(linking='any',descriptor='any')
        chemAtomNmrRef = chemCompVarNmrRef.findFirstChemAtomNmrRef(name=atomName)

        return chemAtomNmrRef.meanValue
