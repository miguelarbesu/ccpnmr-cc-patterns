
def start_cc_diagrams(argServer):
    '''Open the plug-in.'''

    CCPaternsPopup(argServer.parent)

import simulatedPartialSpectrum
reload(simulatedPartialSpectrum)

from simulatedPartialSpectrum import SimulatedPartialSpectrum
from simulatedPartialSpectrum import get_carbons
from itertools import product
from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.popups.BasePopup import BasePopup
from memops.gui.Frame import Frame
from memops.gui.PulldownList import PulldownList
from memops.gui.ButtonList import ButtonList
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.PartitionedSelector import PartitionedSelector


BLUES = ['#a6cee3', '#1f78b4']
GREENS = ['#b2df8a', '#33a02c']
REDS = ['#fb9a99', '#e31a1c']
GREY = '#333333'


class CCPaternsPopup(BasePopup):
    '''The pop-up showing expected peak pattern
       in 13C-13C correlations based on known and
       reference chemical shifts and labeling shemes.
       Also information about assignment status of
       peaks is shown.

    '''

    def __init__(self, parent, *args, **kw):
        '''Init for the pop-up window.'''

        self.font = 'Helvetica 10'
        self.sFont = 'Helvetica %d'
        self.project = parent.project
        self.guiParent = parent
        self.simulated_spectrum = SimulatedPartialSpectrum()
        self.waiting = False
        self.labeling_threshold = 0.1
        self.visible_peaks = []

        BasePopup.__init__(self, parent, title="Carbon Carbon Patterns", **kw)

    def open(self):

        self.updateAfter()
        BasePopup.open(self)

    def body(self, guiFrame):
        '''Setting all the elements of the GUI.'''

        self.geometry('850x750')

        settingsFrame = Frame(guiFrame)
        settingsFrame.grid(row=0, column=0, sticky='nsew')
        residueFrame = Frame(guiFrame)
        residueFrame.grid(row=1, column=0, sticky='nsew')
        diagramFrame = Frame(guiFrame)
        diagramFrame.grid(row=2, column=0, sticky='nsew')

        row = 0
        column = 0

        # Spectrum
        Label(settingsFrame, text='spectrum:', grid=(row, column))
        column += 1
        spectra_pulldown = PulldownList(settingsFrame,
                                        self.update_spectrum,
                                        grid=(row, column))
        spectra = self.get_all_spectra()
        names = [spec.name for spec in spectra]
        spectra_pulldown.setup(names, spectra, 0)
        self.simulated_spectrum.update_spectrum(spectra[0])

        # Peak list
        column += 1
        Label(settingsFrame, text='peak list:', grid=(row, column))
        column += 1
        self.peak_list_pulldown = PulldownList(settingsFrame,
                                               self.simulated_spectrum.update_peak_list,
                                               grid=(row, column))
        self.update_peak_list_pulldown(self.simulated_spectrum.spectrum)


        # Labeling Scheme
        column += 1
        Label(settingsFrame, text='Labeling Scheme:', grid=(row, column))
        column += 1
        scheme_pulldown = PulldownList(settingsFrame,
                                       self.simulated_spectrum.update_labeling_scheme,
                                       grid=(row, column))
        schemes = self.get_all_labeling_schemes()
        names = ['automatic', 'None'] + [scheme.name for scheme in schemes[2:]]
        scheme_pulldown.setup(names, schemes, 0)
        self.simulated_spectrum.update_labeling_scheme('automatic')

        # Shift list
        column += 1
        Label(settingsFrame, text='Shift List:', grid=(row, column))
        column += 1
        shift_list_pulldown = PulldownList(settingsFrame,
                                           self.simulated_spectrum.update_shiftList,
                                           grid=(row, column))
        shift_lists = getShiftLists(self.project.findFirstNmrProject())
        names = ['automatic'] + [shift_list.name for shift_list in shift_lists]
        shift_lists = ['automatic'] + shift_lists
        shift_list_pulldown.setup(names, shift_lists, 0)
        self.simulated_spectrum.update_shiftList('automatic')

        # Molecular Chain
        column += 1
        Label(settingsFrame, text='Chain:', grid=(row, column))
        column += 1
        self.molPulldown = PulldownList(settingsFrame,
                                        self.update_residues,
                                        grid=(row, column))
        self.updateChains()

        row = 0
        column = 0
        texts = ['Previous']
        commands = [self.previous]
        self.previousButton = ButtonList(residueFrame, commands=commands,
                                         texts=texts)
        self.previousButton.grid(row=row, column=column, sticky='ew')

        column += 1
        tipText = '''The Number of the residue in the sequence
                     to display fingerprints for'''
        self.residue_entry1 = IntEntry(residueFrame, grid=(row, column),
                                       width=7, text=1,
                                       returnCallback=self.update_residues,
                                       tipText=tipText)
        column += 1
        Label(residueFrame, text='-', grid=(row, column))
        column += 1
        self.residue_entry2 = IntEntry(residueFrame, grid=(row, column),
                                       width=7, text=2,
                                       returnCallback=self.update_residues,
                                       tipText=tipText)
        column += 1
        texts = ['Next']
        commands = [self.next]
        self.nextButton = ButtonList(residueFrame, commands=commands,
                                     texts=texts)
        self.nextButton.grid(row=row, column=column, sticky='ew')

        column += 1
        texts = ['Show Complete Spectrum']
        commands = [self.select_all_residues]
        show_all_button = ButtonList(residueFrame, commands=commands,
                                     texts=texts)
        show_all_button.grid(row=row, column=column, sticky='ew')

        column += 1
        texts = ['Save peaks']
        commands = [self.save_peaks]
        show_all_button = ButtonList(residueFrame, commands=commands,
                                     texts=texts)
        show_all_button.grid(row=row, column=column, sticky='ew')


        column += 1
        self.peak_text = Label(residueFrame, text=' ', grid=(row, column))

        residueFrame.expandGrid(0, 6)
        column += 1
        self.ppm_text = Label(residueFrame, text=' ', grid=(row, column), sticky='e')

        guiFrame.expandGrid(2, 0)
        diagramFrame.expandGrid(1, 0)
        self.patternSelector = PartitionedSelector(diagramFrame,
                                                   self.toggle_pattern,
                                                   tipText=tipText,
                                                   maxRowObjects=3)
        self.patternSelector.grid(row=0, column=0, columnspan=1, sticky='ew')

        self.diagram = ScrolledGraph(diagramFrame, symbolSize=2, reverseX=True,
                                     reverseY=True, yGrid=False, width=500,
                                     height=500, xLabel='Chemical shift',
                                     yLabel='Chemical shift',
                                     motionCallback=self.update_for_mouse_location)
        self.diagram.grid(row=1, column=0, columnspan=1, sticky='nsew')

        self.simulated_spectrum.register_listener(self.update_diagram)

        self.update_residues()
        #self.select_all_residues()

    def select_all_residues(self):
        chain = self.molPulldown.getObject()
        residues = chain.sortedResidues()
        self.simulated_spectrum.update_residues(residues)
        colors = [REDS[1], BLUES[1], GREENS[1]]

        labels=['all intra-residual', 'all inter-residual', ' ']
        self.patternSelector.update(objects=['intra1', 'seq', 'intra2'],
                                    labels=labels,
                                    colors=colors)

    def previous(self):
        '''Called by previous button. Decreases both
           shown residues by 1.

        '''
        index1 = self.residue_entry1.get()
        index2 = self.residue_entry2.get()
        self.residue_entry1.set(index1 - 1)
        self.residue_entry2.set(index2 - 1)
        self.update_residues()

    def next(self):
        '''Called by next button. Increases both
           shown residues by 1.

        '''
        index1 = self.residue_entry1.get()
        index2 = self.residue_entry2.get()
        self.residue_entry1.set(index1 + 1)
        self.residue_entry2.set(index2 + 1)
        self.update_residues()

    def toggle_pattern(self, residue):
        '''Called when one of the buttons on the
           partitioned selector is pressed toggling
           on and off the sub-patterns corresponding
           either to one of the sets of intra-residual
           peaks or the inter-residual peaks.

        '''
        self.update_diagram()
    
    
    def save_peaks(self):
        '''Save the expected peaks: names, positions, weights.
        '''
        spec = self.simulated_spectrum

        chain = self.molPulldown.getObject()
        chain_name = '%s-%s' % (chain.molSystem.code, chain.code)
        # Automatic scheme not properly defined, drop for now 
        # scheme = spec.labelingScheme
        # fileName = 'expected-peaks_%s_%s.tsv' % (chain_name, scheme)
        fileName = 'expected-peaks_%s.tsv' % (chain_name)
        fileHandle = open(fileName, 'w')
    
        peaks = self.visible_peaks

        for peak in peaks:
            ppm = peak.get_xy()
            weight = peak.colabeling
            fileHandle.write('%.3f\t%.3f\t%.3f\n' % (ppm[0], ppm[1], weight))
        
        print 'Saved expected peaks to %s' % (fileName)
        
        fileHandle.close()
        # The following block each residue name and intra-residue C-C correlations.
        # for residue in spec.residues:
        #     peaks_one_residue = []
        #     carbons = get_carbons(residue)
        #     for carbon1, carbon2 in product(carbons, repeat=2):
        #         print carbon1.getResidue().getCcpCode(), carbon1.getName(), carbon2.getName()

        # fileName = argServer.askString('Output file name', defaultName)
        # if not fileName:
        #     return


    def updatePatternSelector(self):
        '''Change the labels on the partitioned
           selector to match the current residues.

        '''
        colors = [REDS[1], BLUES[1], GREENS[1]]
        res1, res2 = self.get_current_residues()
        labels = ['{} {}'.format(res1.seqCode, res1.ccpCode),
                  'inter-residual',
                  '{} {}'.format(res2.seqCode, res2.ccpCode)]

        self.patternSelector.update(objects=['intra1', 'seq', 'intra2'],
                                    labels=labels,
                                    colors=colors)

    def update_spectrum(self, spectrum):
        '''Update the spectrum that is used.

        '''
        self.update_peak_list_pulldown(spectrum)
        self.simulated_spectrum.update_spectrum(spectrum)

    def update_peak_list_pulldown(self, spectrum):
        '''Update the peak list pulldown. This is necessary
           every time a new spectrum is selected, as they
           have different peak lists.

        '''
        peak_lists = spectrum.sortedPeakLists()
        active = spectrum.getActivePeakList()
        index = peak_lists.index(active)
        texts = ['{}:{}'.format(pl.serial, pl.name) for pl in peak_lists]
        objects = [pl.serial for pl in peak_lists]
        self.peak_list_pulldown.setup(texts, objects, index)
        self.simulated_spectrum.peakListSerial = active.serial

    def update_residues(self, event=None):
        '''Update displayed residues.'''

        residues = self.get_current_residues()
        if residues:
            self.updatePatternSelector()
            self.simulated_spectrum.update_residues(residues)

    def get_current_residues(self):
        '''Get the residues in the chain that
           are selected. Makes sure that the indexes
           are not out of the range of the selected
           chain.

        '''
        residues = []
        chain = self.molPulldown.getObject()
        if not chain or not chain.sortedResidues():
            return residues
        all_residues = chain.sortedResidues()
        for entry in (self.residue_entry1, self.residue_entry2):
            index = entry.get() - 1
            if index < 0:
                index = 0
                entry.set(1)
            elif index >= len(all_residues):
                index = len(all_residues) - 1
                entry.set(len(all_residues))
            residues.append(all_residues[index])
        return residues

    def get_all_labeling_schemes(self):
        '''All labeling schemes appearing in
           the labeling scheme pulldown menu.

        '''
        schemes = ['automatic', None] + self.project.sortedLabelingSchemes()
        return schemes

    def get_all_spectra(self):
        '''All spectra that can be configured.
           Possibly limit this to cc.through-space
           in the future.

        '''
        spectra = []
        nmrProject = self.project.findFirstNmrProject()
        experiments = nmrProject.sortedExperiments()
        for experiment in experiments:
            spectrum = experiment.findFirstDataSource()
            if spectrum:
                spectra.append(spectrum)
        return spectra

    def get_chain_name(self, chain):
        '''Get human readable name for chain.'''

        return '{}:{} ({})'.format(chain.molSystem.code,
                                   chain.code,
                                   chain.molecule.molType)

    def getChains(self):
        '''Get chains from the data model.'''
        chains = []
        if self.project:
            for molSystem in self.project.sortedMolSystems():
                for chain in molSystem.sortedChains():
                    if chain.residues:
                        chains.append(chain)
        return chains

    def updateChains(self, *opt):
        '''Updates the list of chains if a new one is added
           to or deleted from the project. Updates the
           pull down list where a chain can be selected.

        '''
        chains = self.getChains()
        texts = [self.get_chain_name(c) for c in chains]
        self.molPulldown.setup(texts, chains, 0)

    def update_diagram(self):
        '''Redraw the diagram.'''

        selected = self.patternSelector.getSelected()
        spec = self.simulated_spectrum
        all_data_sets = []
        all_symbols = []
        all_colors = []
        self.visible_peaks = []

        intra_peaks = spec.intra_residual_peaks
        if len(spec.intra_residual_peaks) == 2:
            intra1 = spec.intra_residual_peaks[0]
            intra2 = spec.intra_residual_peaks[1]
            inter = spec.sequential_peaks[0]
        elif len(spec.intra_residual_peaks) > 2:
            intra1 = [item for sublist in spec.intra_residual_peaks for item in sublist]
            inter = [item for sublist in spec.sequential_peaks for item in sublist]
            intra2 = []

        for peaks, color_type, kind in [(intra1, REDS, 'intra1'),
                                        (intra2, GREENS, 'intra2'),
                                        (inter, BLUES, 'seq')]:
            if not kind in selected:
                continue
            peaks = self.filter_by_labeling_threshold(peaks)
            data_sets, symbols, colors = self.create_data_sets(peaks,
                                                               color_type)
            all_data_sets.extend(data_sets)
            all_symbols.extend(symbols)
            all_colors.extend(colors)
            self.visible_peaks.extend(peaks)

        # diagonal:
        all_data_sets.append([(-30.0, -30.0, 0.01), (200.0, 200.0, 0.01)])
        all_symbols.append('errorCircle')
        all_colors.append(GREY)
        self.diagram.update(dataSets=all_data_sets,
                            dataColors=all_colors,
                            symbols=all_symbols)
        # print "Here is the data\n", intra_peaks, data_sets

    def create_data_sets(self, peaks, colors):
        '''Make data sets that can be read by the
           scrolled graph. The size of the circles
           is determined by the colabeling of the
           atoms correlated by the peak. The color
           indicates whether a peak with this
           assignment is present in the selected
           peak list.

        '''
        data_sets = []
        symbols = []
        data_colors = []

        for peak in peaks:
            data_point = peak.get_xy() + [peak.colabeling**0.5 * 1.5]
            data_sets.append([data_point])
            symbols.append('errorCircle')
            if peak.assigned_peaks:
                data_colors.append(colors[1])
            else:
                data_colors.append(colors[0])
        return data_sets, symbols, data_colors

    def filter_by_labeling_threshold(self, peaks):
        '''Only return peaks exceeding the
           minimal colabeling.

        '''
        return [peak for peak in peaks if peak.colabeling >= self.labeling_threshold]

    def update_for_mouse_location(self, event):
        '''Update some informational text and
           crosshairs on the spectra.

        '''
        x, y = self.diagram.getPlotCoords(event)
        self.update_crosshairs(x, y, event)
        self.update_peak_text(self.find_closest_peak(x, y))
        self.update_ppm_text(x, y)

    def update_crosshairs(self, x, y, event):
        '''Update the crosshairs on the spectra based
           on the mouse location on the diagram.

        '''
        typeLocation = []
        for panelType in self.analysisProject.panelTypes:
            axisType = panelType.axisType
            if axisType.name == '13C':
                typeLocation.append((panelType, x))
                typeLocation.append((panelType, y))
        self.diagram.mouseEnter(event)
        self.parent.drawWindowCrosshairs(typeLocation)

    def find_closest_peak(self, x, y):
        '''Find the closest peak to a coordinate
           x, y. Not very efficient, but for this
           application fast enough.

        '''
        peaks = self.visible_peaks
        dist_peaks = []
        for peak in peaks:
            peakx, peaky = peak.get_xy()
            dist = ((peakx - x)**2 + (peaky - y)**2)**0.5
            dist_peaks.append((dist, peak))
        if dist_peaks:
            distance, peak = min(dist_peaks)
            if distance < 10:
                return peak

    def update_peak_text(self, peak):
        '''Update an informational text telling
           which to carbons are correlated by
           the crosspeak.

        '''
        if not peak:
            self.peak_text.set('')
            return
        dims = peak.simulatedDims
        data = []
        for dim in dims:
            data.append(dim.atom.residue.seqCode)
            data.append(dim.atom.residue.ccpCode)
            data.append(dim.atom.name)
            if dim.shift_is_assigned:
                data.append('assigned')
            else:
                data.append('unassigned')

        text = '{}{} {} {} - {}{} {} {}'.format(*data)
        self.peak_text.set(text)

    def update_ppm_text(self, x, y):
        '''Showing the current x,y coordinates.'''

        #self.ppm_text.set('{0:.2f}, {0:.2f}'.format(x, y))
        self.ppm_text.set('{0:.2f}, {1:.2f}'.format(x, y))
