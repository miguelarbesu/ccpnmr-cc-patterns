
import simulatedPeakDimension
reload(simulatedPeakDimension)
from simulatedPeakDimension import SimulatedPeakDimension
from ccp.util.LabeledMolecule import getExperimentAtomPairFractions, getSchemeAtomPairFractions


class SimulatedPeak(object):
    """Peak that is expected to appear in a spectrum."""
    def __init__(self, simulated_partial_spectrum=None, atoms=None):
        self.simulated_partial_spectrum = simulated_partial_spectrum
        self.simulatedDims = [SimulatedPeakDimension(simulatedPeak=self, atom=atom) for atom in atoms]
        self.colabeling = 1.0
        self.update_colabeling()
        self.assigned_peaks = set()
        self.update_assigned_peaks()

    def get_xy(self):
        '''In case of 2D spectrum returns x,y
           coordinate in ppm of the peak.
           In nD spectra list with n coordinates.

        '''
        return [dim.shift for dim in self.simulatedDims]

    def update_position(self):
        '''Recalculate the expected position
           in all dimensions of the peak.

        '''
        for dim in self.simulatedDims:
            dim.update_shift()

    def update_assigned_peaks(self):
        '''Updates list of peaks that are assigned_peaks
           in a peak list to same correlation as this
           expected peak.

        '''
        self.assigned_peaks = self.get_assigned_peaks()

    def get_assigned_peaks(self):
        '''Get all peaks in a peak list that are
           assigned to the same correlation between
           atoms as represented by this expected peak.

        '''
        atoms = [dim.atom for dim in self.simulatedDims]
        spectrum = self.simulated_partial_spectrum.spectrum
        if not spectrum:
            return set()
        peakListSerial = self.simulated_partial_spectrum.peakListSerial
        contribs1 = get_peakDimContribs(atoms[0], dim_number=1,
                                        spectrum=spectrum,
                                        peakListSerial=peakListSerial)
        contribs2 = get_peakDimContribs(atoms[1], dim_number=2,
                                        spectrum=spectrum,
                                        peakListSerial=peakListSerial)
        peaks1 = set([contrib.peakDim.peak for contrib in contribs1])
        peaks2 = set([contrib.peakDim.peak for contrib in contribs2])
        peaks = peaks1.intersection(peaks2)
        return peaks

    def update_colabeling(self):
        '''Update the colabeling of the two atoms
           that are correlated by this peak. Only
           works correctly for 2D.

        '''
        atom1 = self.simulatedDims[0].atom
        atom2 = self.simulatedDims[1].atom

        if self.simulated_partial_spectrum.labelingScheme is None:
            self.colabeling = 1.0
            return
        elif self.simulated_partial_spectrum.labelingScheme == 'automatic':
            experiment = self.simulated_partial_spectrum.spectrum.experiment
            fracs = getExperimentAtomPairFractions(experiment, atom1, atom2)
        else:
            scheme = self.simulated_partial_spectrum.labelingScheme
            fracs = getSchemeAtomPairFractions(scheme, atom1, atom2)
        if not fracs:
            self.colabeling = 1.0
        else:
            fraction = fracs[('13C', '13C')]
            self.colabeling = fraction


def get_peakDimContribs(atom, dim_number=None, spectrum=None,
                        peakListSerial=None):
    '''Get all contributions to dimensions of peaks for
       an atom in a peak list.

    '''

    dimContribs = set()
    if atom.atomSet and atom.atomSet.resonanceSets:
        for resonanceSet in atom.atomSet.resonanceSets:
            resonances = resonanceSet.resonances
            for resonance in resonances:
                if dim_number:
                    dimContribs |= resonance.findAllPeakDimContribs(dim=dim_number)
                else:
                    dimContribs |= resonance.findAllPeakDimContribs()
    if spectrum:
        dimContribs = set([contrib for contrib in dimContribs if contrib.peakDim.peak.peakList.dataSource is spectrum])
    if peakListSerial:
        dimContribs = set([contrib for contrib in dimContribs if contrib.peakDim.peak.peakList.serial == peakListSerial])
    return dimContribs

