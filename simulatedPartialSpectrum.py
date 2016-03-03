
from simulatedPeak import SimulatedPeak
from itertools import product
import simulatedPeak
reload(simulatedPeak)


def test_simulating_spectrum(argServer):
    '''Quick test without GUI.'''

    project = argServer.getProject()
    nmrProject = project.findFirstNmrProject()
    chain = project.findFirstMolSystem(name='MS1').findFirstChain()
    residues = chain.sortedResidues()[199:201]
    experiment = nmrProject.findFirstExperiment(name='029_PDSD400_13C')
    spectrum = experiment.findFirstDataSource()
    sim = SimulatedPartialSpectrum(residues=residues, spectrum=spectrum,
                                   shiftList='automatic',
                                   labelingScheme='automatic',
                                   peakListSerial=1)

    for peak in sim.get_all_peaks():
        print 'position', peak.get_xy()
        print 'labeling', peak.colabeling
        print 'assigned', peak.assigned_peaks


class SimulatedPartialSpectrum(object):
    """Holding all information about peaks in a
       sub-pattern concerning only a few residues.

    """
    def __init__(self, residues=None, labelingScheme=None,
                 spectrum=None, peakListSerial=None, shiftList=None):

        self.residues = residues
        self.labelingScheme = labelingScheme
        self.spectrum = spectrum
        self.peakListSerial = peakListSerial
        self.intra_residual_peaks = [[], []]
        self.sequential_peaks = [[]]
        self.shiftList = shiftList
        self.update_shiftList(shiftList)
        self.listeners = []
        self.predict_peaks()

    def register_listener(self, listener):
        '''Register an observer function to be called
           when peak pattern changes.

        '''
        self.listeners.append(listener)

    def notify_listeners(self):
        '''Call all observer functions.'''
        for listener in self.listeners:
            listener()

    def predict_peaks(self):
        '''Predict intra- and inter-residual
           cross-peak pattern.

        '''
        if not self.residues:
            return
        self.predict_intra_peaks()
        self.predict_sequential_peaks()
        self.notify_listeners()

    def get_all_peaks(self):
        '''Get a flat list with all the peaks.'''

        all_peak_groups = self.intra_residual_peaks + self.sequential_peaks
        for peak_group in all_peak_groups:
            for peak in peak_group:
                yield peak

    def predict_intra_peaks(self):
        '''Predict intra-residual peaks'''
        self.intra_residual_peaks = []
        for residue in self.residues:
            peaks_one_residue = []
            carbons = get_carbons(residue)
            for carbon1, carbon2 in product(carbons, repeat=2):
                peak = SimulatedPeak(simulated_partial_spectrum=self,
                                     atoms=(carbon1, carbon2))
                peaks_one_residue.append(peak)
            self.intra_residual_peaks.append(peaks_one_residue)

    def predict_sequential_peaks(self):
        '''Predict inter-residual peaks'''
        self.sequential_peaks = []
        for res1, res2 in zip(self.residues[:-1], self.residues[1:]):
            peaks_one_residue_pair = []
            carbons1 = get_carbons(res1)
            carbons2 = get_carbons(res2)
            for carbon1, carbon2 in product(carbons1, carbons2):
                peak = SimulatedPeak(simulated_partial_spectrum=self,
                                     atoms=(carbon1, carbon2))
                peaks_one_residue_pair.append(peak)
                peak = SimulatedPeak(simulated_partial_spectrum=self,
                                     atoms=(carbon2, carbon1))
                peaks_one_residue_pair.append(peak)
            self.sequential_peaks.append(peaks_one_residue_pair)

    def update_residues(self, residues):
        '''Update residues and re-predict peak pattern.'''
        self.residues = residues
        self.predict_peaks()

    def update_spectrum(self, spectrum):
        '''Update spectrum and update all
           information on the peaks.

        '''
        if spectrum is self.spectrum:
            return
        self.spectrum = spectrum
        if self.shiftList == 'automatic':
            self.update_shiftList('automatic', no_notify=True)
        if self.labelingScheme == 'automatic':
            self.update_labeling_scheme('automatic', no_notify=True)
        for peak in self.get_all_peaks():
            peak.update_assigned_peaks()
        self.notify_listeners()

    def update_peak_list(self, serial, no_notify=False):
        '''Update peak list.'''
        self.peakListSerial = serial
        for peak in self.get_all_peaks():
            peak.update_assigned_peaks()
        if not no_notify:
            self.notify_listeners()

    def update_labeling_scheme(self, labelingScheme, no_notify=False):
        '''Update labeling scheme.'''
        self.labelingScheme = labelingScheme
        for peak in self.get_all_peaks():
            peak.update_colabeling()
        if not no_notify:
            self.notify_listeners()

    def update_shiftList(self, shiftList, no_notify=False):
        '''Update shift list.'''
        if shiftList == 'automatic':
            new_shiftList = self.spectrum.experiment.shiftList
        else:
            new_shiftList = shiftList
        if self.shiftList is new_shiftList:
            return

        self.shiftList = new_shiftList
        for peak in self.get_all_peaks():
            peak.update_position()
        if not no_notify:
            self.notify_listeners()


def get_carbons(residue):
    '''get all carbon atoms within a residue.'''
    carbons = []
    for atom in residue.sortedAtoms():
        if atom.chemAtom.elementSymbol == 'C':
            carbons.append(atom)
    return carbons
