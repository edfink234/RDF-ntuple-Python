from MakeRDF import MakeRDF
from time import time
import ROOT

#ROOT.gSystem.Load("./MakeRDF_cxx.so")

"""
Declaring the functions that we use for Define and Filter
in this analysis. The bool functions are used in calls to
Filter, and the rest are used in calls to Define.
"""
ROOT.gInterpreter.Declare("""
bool photon_selection(Photon& photon)
{
    if (!photon.photon_id)
    {
        return false;
    }
    
    if (photon.photon_pt < 2500)
    {
        return false;
    }
    
    if (abs(photon.photon_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(photon.photon_eta)) && (abs(photon.photon_eta) < 1.52))
    {
        return false;
    }
    return true;
}

bool truth_photon_selection(TruthParticle& truth_photon)
{
    if (truth_photon.mc_pt < 2500)
    {
        return false;
    }
    
    if (abs(truth_photon.mc_eta) > 2.37)
    {
        return false;
    }
    
    if ((1.37 < abs(truth_photon.mc_eta)) && (abs(truth_photon.mc_eta) < 1.52))
    {
        return false;
    }
    return true;
}

bool track_selection(Track& track)
{
    if (track.track_pt < 100)
    {
        return false;
    }
    
    if (abs(track.track_eta) > 2.5)
    {
        return false;
    }
    
    return true;
}

RVec<Photon> selected_photons_func(RVec<Photon> photons)
{
    photons.erase(std::remove_if(photons.begin(),photons.end(),
    [](Photon& x)
    {
        return !photon_selection(x);
        
    }), photons.end());
    
    return photons;
}

RVec<Track> selected_tracks_func(RVec<Track> tracks)
{
   tracks.erase(std::remove_if(tracks.begin(),tracks.end(),
   [](Track& x)
   {
       return !track_selection(x);
       
   }), tracks.end());
   
   return tracks;
}

RVec<TruthParticle> truth_photons_func(RVec<TruthParticle> truth_particles)
{
    truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
    [](TruthParticle& x)
    {
        return (x.mc_pdg_id!=22 || !truth_photon_selection(x));
        
    }), truth_particles.end());
    
    return truth_particles;
}

bool truth_candidates_func(RVec<TruthParticle>& truth_photons)
{
    if (truth_photons.size() == 2)
    {
        PtEtaPhiEVector four_momentum =
        truth_photons[0].Vector() + truth_photons[1].Vector();
        
        if (four_momentum.M() > 5000)
        {
            return true;
        }
    }
    return false;
}

RVec<float> truth_candidates_pt_func(RVec<TruthParticle>& truth_photons)
{
    return RVec<float>({static_cast<float>(truth_photons[0].mc_pt/1e3), static_cast<float>(truth_photons[1].mc_pt/1e3)});
}

ROOT::Math::PtEtaPhiEVector::Scalar truth_candidates_mass_func(RVec<TruthParticle>& diphotons)
{
    PtEtaPhiEVector four_momentum =
    diphotons[0].Vector() + diphotons[1].Vector();
    
    return four_momentum.M()/1e3;
}

RVec<float> diphotons_pt_func(RVec<Photon>& diphotons)
{
    return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
}

ROOT::Math::PtEtaPhiEVector::Scalar diphotons_mass_func(RVec<Photon>& diphotons)
{
    PtEtaPhiEVector four_momentum =
    diphotons[0].Vector() + diphotons[1].Vector();
    
    return four_momentum.M()/1e3;
}

RVec<float> selected_tracks_pt_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_pt;
    selected_tracks_pt.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_pt.push_back(track.track_pt/1e3);
    }
    return selected_tracks_pt;
}

RVec<float> selected_tracks_eta_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_eta;
    selected_tracks_eta.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_eta.push_back(track.track_eta);
    }
    return selected_tracks_eta;
}

bool no_tracks_inv_mass_func(RVec<Track>& selected_tracks, RVec<Photon>& diphotons)
{
    PtEtaPhiEVector four_momentum =
    diphotons[0].Vector() + diphotons[1].Vector();
    
    return (selected_tracks.empty() && four_momentum.M() > 5000);
}

RVec<float> diphotons_pt_no_tracks_inv_mass_func(RVec<Photon>& diphotons)
{
    return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
}

ROOT::Math::PtEtaPhiEVector::Scalar diphotons_mass_no_tracks_inv_mass_func(RVec<Photon>& diphotons)
{
    PtEtaPhiEVector four_momentum =
    diphotons[0].Vector() + diphotons[1].Vector();
    
    return four_momentum.M()/1e3;
}

RVec<float> selected_tracks_pt_no_tracks_inv_mass_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_pt;
    selected_tracks_pt.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_pt.push_back(track.track_pt/1e3);
    }
    return selected_tracks_pt;
}

RVec<float> selected_tracks_eta_no_tracks_inv_mass_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_eta;
    selected_tracks_eta.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_eta.push_back(track.track_eta);
    }
    return selected_tracks_eta;
    
}

RVec<float> diphotons_pt_no_tracks_func(RVec<Photon>& diphotons)
{
    return RVec<float>({static_cast<float>(diphotons[0].photon_pt/1e3), static_cast<float>(diphotons[1].photon_pt/1e3)});
}

ROOT::Math::PtEtaPhiEVector::Scalar diphotons_mass_no_tracks_func(RVec<Photon>& diphotons)
{
    PtEtaPhiEVector four_momentum =
    diphotons[0].Vector() + diphotons[1].Vector();
    
    return four_momentum.M()/1e3;
}

RVec<float> selected_tracks_pt_no_tracks_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_pt;
    selected_tracks_pt.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_pt.push_back(track.track_pt/1e3);
    }
    return selected_tracks_pt;
    
}

RVec<float> selected_tracks_eta_no_tracks_func(RVec<Track>& tracks)
{
    RVec<float> selected_tracks_eta;
    selected_tracks_eta.reserve(tracks.size());
    for (auto& track: tracks)
    {
        selected_tracks_eta.push_back(track.track_eta);
    }
    return selected_tracks_eta;
    
}

""")

if __name__ == "__main__":
    start = time()
    files = ("/Users/edwardfinkelstein/ATLAS_axion/user.kschmied.28655874._000025.LGNTuple.root",)
    ROOT.Event.systematics = {"EG_RESOLUTION_ALL__1down"};
    df = MakeRDF(files) #Create the dataframe
#    df.Describe().Print()
    
#    An example of filter. In this case, no events are filtered
#    since mc is true
    firstTriggerCut = df.Define("mc", "true").Filter("mc || std::find(trigger_passed_triggers.begin(), trigger_passed_triggers.end(), \"HLT_hi_upc_FgapAC3_hi_gg_upc_L1TAU1_TE4_VTE200\") != trigger_passed_triggers.end()")
    
#    print(firstTriggerCut.Count().GetValue())

#    Book columns for selected photons and selected tracks using
#    functions defined above
    selected_photons_and_tracks = firstTriggerCut.Define("selected_photons","selected_photons_func(photons)")\
                       .Define("selected_tracks", "selected_tracks_func(tracks)")
    
#    Books columns for selected truth_particles
    truth_photons = firstTriggerCut.Define("truth_photons", "truth_photons_func(truth_particles)")
    
#    truth_photons.Describe().Print()
#    Book a filtration of truth_photons; for the events that pass
#    book columns for the p_t and mass
    truth_candidates = truth_photons.Filter("truth_candidates_func(truth_photons)")\
                                    .Define("truth_candidates_pt", "truth_candidates_pt_func(truth_photons)")\
                                    .Define("truth_candidates_mass", "truth_candidates_mass_func(truth_photons)")
    
#    Filter events that don't have exactly two photons; for the events
#    that pass, define columns for the pt, mass, number of selected tracks
#    and the pt and eta of the selected tracks
    diphotons_and_tracks = selected_photons_and_tracks.Filter("selected_photons.size() == 2").Define("diphotons_pt", "diphotons_pt_func(selected_photons)")\
           .Define("diphotons_mass", "diphotons_mass_func(selected_photons)")\
           .Define("num_tracks", "selected_tracks.size()")\
           .Define("selected_tracks_pt", "selected_tracks_pt_func(selected_tracks)")\
           .Define("selected_tracks_eta", "selected_tracks_eta_func(selected_tracks)")
    
#    Filter events that have tracks or diphoton inv mass <= 5 GeV,
#    then define columns for the diphoton pt, inv mass, number of
#    tracks, selected tracks pt and eta.
    no_tracks_inv_mass = diphotons_and_tracks.Filter("no_tracks_inv_mass_func(selected_tracks, selected_photons)")\
           .Define("diphotons_pt_no_tracks_inv_mass", "diphotons_pt_no_tracks_inv_mass_func(selected_photons)")\
           .Define("diphotons_mass_no_tracks_inv_mass", "diphotons_mass_no_tracks_inv_mass_func(selected_photons)")\
           .Define("num_tracks_no_tracks_inv_mass", "selected_tracks.size()")\
           .Define("selected_tracks_pt_no_tracks_inv_mass", "selected_tracks_pt_no_tracks_inv_mass_func(selected_tracks)")\
           .Define("selected_tracks_eta_no_tracks_inv_mass", "selected_tracks_eta_no_tracks_inv_mass_func(selected_tracks)")
         
#    Filter events that have tracks, then define columns for diphoton pt,
#    inv mass, number of tracks, selected tracks pt and eta.
    no_tracks = diphotons_and_tracks.Filter("selected_tracks.empty()")\
           .Define("diphotons_pt_no_tracks","diphotons_pt_no_tracks_func(selected_photons)")\
           .Define("diphotons_mass_no_tracks","diphotons_mass_no_tracks_func(selected_photons)")\
           .Define("num_tracks_no_tracks","selected_tracks.size()")\
           .Define("selected_tracks_pt_no_tracks","selected_tracks_pt_no_tracks_func(selected_tracks)")\
           .Define("selected_tracks_eta_no_tracks","selected_tracks_eta_no_tracks_func(selected_tracks)")
    
#    A list of all of the Histograms we want to create. The event loop still
#    hasn't been triggered yet!
    histos = (
        truth_candidates.Histo1D(("TruthRecoPhotonPt", "TruthRecoPhotonPt", 20, 0, 25), "truth_candidates_pt"),
        truth_candidates.Histo1D(("TruthRecoPhotonPtTight", "TruthRecoPhotonPtTight", 100, 0, 10), "truth_candidates_pt"),
        truth_candidates.Histo1D(("TruthCandidateMass", "TruthCandidateMass", 250, 0, 50), "truth_candidates_mass"),
        truth_candidates.Histo1D(("TruthCandidateMassLarge", "TruthCandidateMassLarge", 200, 0, 200), "truth_candidates_mass"),
        truth_candidates.Histo1D(("TruthCandidateMassLargeFine", "TruthCandidateMassLargeFine", 600, 0, 200), "truth_candidates_mass"),
        
        diphotons_and_tracks.Histo1D(("00NoCutsRecoPhotonPt", "00NoCutsRecoPhotonPt", 20, 0, 25), "diphotons_pt"),
        diphotons_and_tracks.Histo1D(("00NoCutsRecoPhotonPtTight", "00NoCutsRecoPhotonPtTight", 100, 0, 10), "diphotons_pt"),
        diphotons_and_tracks.Histo1D(("00NoCutsCandidateMass", "00NoCutsCandidateMass", 250, 0, 50), "diphotons_mass"),
        diphotons_and_tracks.Histo1D(("00NoCutsCandidateMassLarge", "00NoCutsCandidateMassLarge", 200, 0, 200), "diphotons_mass"),
        diphotons_and_tracks.Histo1D(("00NoCutsCandidateMassLargeFine", "00NoCutsCandidateMassLargeFine", 600, 0, 200), "diphotons_mass"),
        diphotons_and_tracks.Histo1D(("00NoCutsTrackingNumTracks", "00NoCutsTrackingNumTracks", 20, 0, 20), "num_tracks"),
        diphotons_and_tracks.Histo1D(("00NoCutsTrackingTrackPt", "00NoCutsTrackingTrackPt", 140, 0, 7), "selected_tracks_pt"),
        diphotons_and_tracks.Histo1D(("00NoCutsTrackingTrackEta", "00NoCutsTrackingTrackEta", 50, -2.5, 2.5), "selected_tracks_eta"),

        no_tracks_inv_mass.Histo1D(("02MassCutRecoPhotonPt", "02MassCutRecoPhotonPt", 20, 0, 25), "diphotons_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutRecoPhotonPtTight", "02MassCutRecoPhotonPtTight", 100, 0, 10), "diphotons_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutCandidateMass", "02MassCutCandidateMass", 250, 0, 50), "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutCandidateMassLarge", "02MassCutCandidateMassLarge", 200, 0, 200), "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutCandidateMassLargeFine", "02MassCutCandidateMassLargeFine", 600, 0, 200), "diphotons_mass_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutTrackingNumTracks", "02MassCutTrackingNumTracks", 20, 0, 20), "num_tracks_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutTrackingTrackPt", "02MassCutTrackingTrackPt", 140, 0, 7), "selected_tracks_pt_no_tracks_inv_mass"),
        no_tracks_inv_mass.Histo1D(("02MassCutTrackingTrackEta", "02MassCutTrackingTrackEta", 50, -2.5, 2.5), "selected_tracks_eta_no_tracks_inv_mass"),

        no_tracks.Histo1D(("01NoTracksRecoPhotonPt", "01NoTracksRecoPhotonPt", 20, 0, 25), "diphotons_pt_no_tracks"),
        no_tracks.Histo1D(("01NoTracksRecoPhotonPtTight", "01NoTracksRecoPhotonPtTight", 100, 0, 10), "diphotons_pt_no_tracks"),
        no_tracks.Histo1D(("01NoTracksCandidateMass", "01NoTracksCandidateMass", 250, 0, 50), "diphotons_mass_no_tracks"),
        no_tracks.Histo1D(("01NoTracksCandidateMassLarge", "01NoTracksCandidateMassLarge", 200, 0, 200), "diphotons_mass_no_tracks"),
        no_tracks.Histo1D(("01NoTracksCandidateMassLargeFine", "01NoTracksCandidateMassLargeFine", 600, 0, 200), "diphotons_mass_no_tracks"),
        no_tracks.Histo1D(("01NoTracksTrackingNumTracks", "01NoTracksTrackingNumTracks", 20, 0, 20), "num_tracks_no_tracks"),
        no_tracks.Histo1D(("01NoTracksTrackingTrackPt", "01NoTracksTrackingTrackPt", 140, 0, 7), "selected_tracks_pt_no_tracks"),
        no_tracks.Histo1D(("01NoTracksTrackingTrackEta", "01NoTracksTrackingTrackEta", 50, -2.5, 2.5), "selected_tracks_eta_no_tracks")
    )
    
    resultmaps = []
    histNames = []
    for h in histos:
#        Appending a map that contains the nominal and varied Histograms
#        for the systematics in Event.systematics to the list resultmaps
        resultmaps.append(ROOT.RDF.Experimental.VariationsFor(h))
        histNames.append(h.GetName())
    
    for i in range(len(resultmaps)): #For map of histograms
        for var in resultmaps[i].GetKeys(): #For each varied histogram in the map
            c1 = ROOT.TCanvas("","",800, 700) #Create a new canvas
            resultmaps[i][var].Draw("same") #Draw the histogram
            string = str(var).replace(":","")
            
            c1.SaveAs((string+histNames[i]+".png")) #Save the histogram.
            if histNames[i] == "02MassCutTrackingNumTracks":
                print(f"# of events for {var} = {resultmaps[i][var].GetEntries()}")
    
    end = time() #done
    print(f"Time taken = {end-start:0.3f} seconds")
    
    
        

    
