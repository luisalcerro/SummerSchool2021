#!/usr/bin/env python3

"""
  Class to analyze a single JETSCAPE output file

  Author: James Mulligan (james.mulligan@berkeley.edu)
  """

from __future__ import print_function

# General
import sys
import os
import argparse
import yaml
import pyhepmc_ng
import math
# Fastjet via python (from external library heppy)
import fastjet as fj
import ROOT

sys.path.append('../..')
from jetscape_analysis.analysis import analyze_events_base

################################################################
class AnalyzeJetscapeEvents_Example(analyze_events_base.AnalyzeJetscapeEvents_Base):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, config_file='', input_file='', output_dir='', **kwargs):
        super(AnalyzeJetscapeEvents_Example, self).__init__(config_file=config_file,
                                                            input_file=input_file,
                                                            output_dir=output_dir,
                                                            **kwargs)
        self.initialize_user_config()
        print(self)

    # ---------------------------------------------------------------
    # Initialize config file into class members
    # ---------------------------------------------------------------
    def initialize_user_config(self):

        # Read config file
        with open(self.config_file, 'r') as stream:
            config = yaml.safe_load(stream)
        
        self.min_track_pt = config['min_track_pt']
        self.jetR_list = config['jetR']
        self.min_jet_pt = config['min_jet_pt']
        self.n_event_max = config['n_event_max']

    # ---------------------------------------------------------------
    # Initialize output objects
    # ---------------------------------------------------------------
    def initialize_user_output_objects(self):

        # Hadron histograms
        hname = 'hChHadronPt'
        h = ROOT.TH1F(hname, hname, 100, 0, 100)
        h.Sumw2()
        setattr(self, hname, h)
        

        # Lepton histograms

        hname = 'hNLep'
        h = ROOT.TH1F(hname, hname, 10, 0, 10)
        setattr(self, hname, h)

                        
        hname = 'hLepPt'
        h = ROOT.TH1F(hname, hname, 20, 0, 200)
        setattr(self, hname, h)

        hname = 'hLepEta'
        h = ROOT.TH1F(hname, hname, 50, -10, 10)
        setattr(self, hname, h)

        hname = 'hLepPhi'
        h = ROOT.TH1F(hname, hname, 50, -4, 4)
        setattr(self, hname, h)

        hname = 'hW2massLep'
        h = ROOT.TH1F(hname, hname, 50, 60, 100)
        setattr(self, hname, h)

       # Jet histograms
        for jetR in self.jetR_list:

            hname = 'hJetPt_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, 8, 80, 180)
            setattr(self, hname, h)

            hname = 'hJetEta_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, 50, -10, 10)
            setattr(self, hname, h)

            hname = 'hJetPhi_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, 50, -1, 7)
            setattr(self, hname, h)
                                    
            hname = 'hNJet_R{}'.format(jetR)
            h = ROOT.TH1F(hname, hname, 10, 0, 10)
            setattr(self, hname, h)

            

    # ---------------------------------------------------------------
    # Analyze a single event -- fill user-defined output objects
    # ---------------------------------------------------------------
    def analyze_event(self, event):

        # Get list of hadrons from the event, and fill some histograms
        hadrons = event.hadrons(min_track_pt=self.min_track_pt)
        self.fill_hadron_histograms(hadrons)
	
	    # Get list of leptons from the event, and fill some histograms
        leptons = event.leptons(min_track_pt=self.min_track_pt)
        self.fill_lepton_histograms(leptons)
        
        # Create list of fastjet::PseudoJets
        fj_hadrons = self.fill_fastjet_constituents(hadrons)

        # Loop through specified jet R
        for jetR in self.jetR_list:

            # Set jet definition and a jet selector
            jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
            jet_selector = fj.SelectorPtMin(self.min_jet_pt) & fj.SelectorAbsRapMax(10.)
            if self.debug_level > 0:
                print('jet definition is:', jet_def)
                print('jet selector is:', jet_selector, '\n')

            # Do jet finding
            cs = fj.ClusterSequence(fj_hadrons, jet_def)
            jets = fj.sorted_by_pt(cs.inclusive_jets())
            jets_selected = jet_selector(jets)

            # Fill some jet histograms
            self.fill_jet_histograms(jets_selected, jetR)

    # ---------------------------------------------------------------
    # Fill hadron histograms
    # ---------------------------------------------------------------
    def fill_hadron_histograms(self, hadrons):
        # Loop through hadrons
        for hadron in hadrons:
            
            # Fill some basic hadron info
            pid = hadron.pid
            pt = hadron.momentum.pt()
            eta = hadron.momentum.eta()
            phi = hadron.momentum.phi()
	    
            # Fill charged hadron histograms (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
            if abs(pid) in [211, 321, 2212, 3222, 3112, 3312, 3334]:
                getattr(self, 'hChHadronPt').Fill(pt, 1/pt) # Fill with weight 1/pt, to form 1/pt dN/dpt

    
    # ---------------------------------------------------------------
    # Fill lepton histograms
    # ---------------------------------------------------------------
    def fill_lepton_histograms(self, leptons):
        # Loop through leptons
        Nlep=0
        muptMax = 20.0
        mubarptMax = 20.0
        numubarptMax = 0.0
        numuptMax = 0.0

        muMax = ''
        mubarMax = ''
        numuMax = ''
        numubarMax = ''
        for lepton in leptons:
            # Fill some basic hadron info
            pid = lepton.pid
            pt = lepton.momentum.pt()
            eta = lepton.momentum.eta()
            phi = lepton.momentum.phi()

            # Fill leptons    
            if abs(pid) == 13:
                Nlep +=1
                getattr(self, 'hLepPt').Fill(pt, 1.0/self.n_event_max)
                getattr(self, 'hLepEta').Fill(eta)
                getattr(self, 'hLepPhi').Fill(phi)

            # W mass
            if pid == 13:
                if pt>muptMax:
                    muptMax = pt
                    muMax = lepton.momentum
            if pid == -14:
                if pt>numubarptMax:
                    munubarptMax = pt
                    numubarMax = lepton.momentum

            if pid == -13:
                if pt>mubarptMax:
                    mubarptMax = pt
                    mubarMax = lepton.momentum
            if pid == 14:
                if pt>numuptMax:
                    numuptMax = pt
                    numuMax = lepton.momentum

        getattr(self, 'hNLep').Fill(Nlep)

        W2massLep1 = ''
        W2massLep2 = ''
        if(bool(muMax) and bool(numubarMax)):
            W2massLep1 = (muMax+numubarMax).m()
        if(bool(mubarMax) and bool(numuMax)):
            W2massLep2 = (mubarMax+numuMax).m()
        if(bool(W2massLep1)):getattr(self, 'hW2massLep').Fill(W2massLep1)
        if(bool(W2massLep2)): getattr(self, 'hW2massLep').Fill(W2massLep2)
    # ---------------------------------------------------------------
    # Fill jet histograms
    # ---------------------------------------------------------------
    def fill_jet_histograms(self, jets, jetR):

        for jet in jets:
            jet_pt = jet.pt()
            jet_eta = jet.eta()
            jet_phi = jet.phi()
            getattr(self, 'hJetPt_R{}'.format(jetR)).Fill(jet_pt, 1.0/self.n_event_max)
            getattr(self, 'hJetEta_R{}'.format(jetR)).Fill(jet_eta)
            getattr(self, 'hJetPhi_R{}'.format(jetR)).Fill(jet_phi)
        getattr(self, 'hNJet_R{}'.format(jetR)).Fill(len(jets))

##################################################################
if __name__ == "__main__":
    # Define arguments
    parser = argparse.ArgumentParser(description="Generate JETSCAPE events")
    parser.add_argument(
        "-c",
        "--configFile",
        action="store",
        type=str,
        metavar="configFile",
        default="/home/jetscape-user/JETSCAPE-analysis/config/jetscapeAnalysisConfig.yaml",
        help="Path of config file for analysis",
    )
    parser.add_argument(
        "-i",
        "--inputDir",
        action="store",
        type=str,
        metavar="inputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Input directory containing JETSCAPE output files",
    )
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        type=str,
        metavar="outputDir",
        default="/home/jetscape-user/JETSCAPE-analysis/TestOutput",
        help="Output directory for output to be written to",
    )

    # Parse the arguments
    args = parser.parse_args()

    # If invalid configFile is given, exit
    if not os.path.exists(args.configFile):
        print('File "{0}" does not exist! Exiting!'.format(args.configFile))
        sys.exit(0)

    # If invalid inputDir is given, exit
    if not os.path.exists(args.inputDir):
        print('File "{0}" does not exist! Exiting!'.format(args.inputDir))
        sys.exit(0)

    analysis = AnalyzeJetscapeEvents_Example(config_file=args.configFile, input_dir=args.inputDir, output_dir=args.outputDir)
    analysis.analyze_jetscape_events()
