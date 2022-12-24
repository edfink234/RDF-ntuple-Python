from MakeRDF import MakeRDF
from time import time
import ROOT

#ROOT.gSystem.Load("./MakeRDF_cxx.so")

ROOT.gInterpreter.Declare("""

RVec<TruthParticle> stable_truth_photons_func(RVec<TruthParticle> truth_particles)
{
    truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
    [](const TruthParticle& x)
    {
        return (x.mc_pdg_id!=22 || x.mc_status!=1);
        
    }), truth_particles.end());
    
    return truth_particles;
}

RVec<float> stable_truth_photons_pt_func(const RVec<TruthParticle>& stable_truth_photons)
{
    RVec<float> pt;
    pt.reserve(stable_truth_photons.size());
    for (auto &i: stable_truth_photons)
    {
        pt.push_back(i.mc_pt/1e3);
    }
    return pt;
}

RVec<float> stable_truth_photons_eta_func(const RVec<TruthParticle>& stable_truth_photons)
{
    RVec<float> eta;
    eta.reserve(stable_truth_photons.size());
    for (auto &i: stable_truth_photons)
    {
        eta.push_back(i.mc_eta);
    }
    return eta;
}

RVec<TruthParticle> stable_truth_leptons_func(RVec<TruthParticle> truth_particles)
{
    truth_particles.erase(std::remove_if(truth_particles.begin(),truth_particles.end(),
    [](const TruthParticle& x)
    {
          return (((abs(x.mc_pdg_id)!=11)
              && (abs(x.mc_pdg_id)!=12)
              && (abs(x.mc_pdg_id)!=13)
              && (abs(x.mc_pdg_id)!=14)
              && (abs(x.mc_pdg_id)!=15)
              && (abs(x.mc_pdg_id)!=16)
              && (abs(x.mc_pdg_id)!=17)
              && (abs(x.mc_pdg_id)!=18))
              || (x.mc_status!=1));

    }), truth_particles.end());
    
    return truth_particles;
}

using namespace ROOT::Math::VectorUtil; // DeltaR
bool passed_truth_leptons_func(const RVec<TruthParticle>& stable_truth_leptons)
{
    return
    (
        (stable_truth_leptons.size()==2)
     
        &&
     
        (DeltaR(
            PtEtaPhiEVector(stable_truth_leptons[0].mc_pt,
                            stable_truth_leptons[0].mc_eta,
                            stable_truth_leptons[0].mc_phi,
                            stable_truth_leptons[0].mc_e),

            PtEtaPhiEVector(stable_truth_leptons[1].mc_pt,
                            stable_truth_leptons[1].mc_eta,
                            stable_truth_leptons[1].mc_phi,
                            stable_truth_leptons[1].mc_e)
           ) > 0.01)
     
        &&
     
        (stable_truth_leptons[0].mc_charge == -1*stable_truth_leptons[1].mc_charge)

        &&

        ((((stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e)
        * (stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e))
        -
        ((stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)
        * (stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)))
        >= 6561e6)

        &&

        ((((stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e)
        * (stable_truth_leptons[0].mc_e + stable_truth_leptons[1].mc_e))
        -
        ((stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)
        * (stable_truth_leptons[0].mc_pt + stable_truth_leptons[1].mc_pt)))
        <= 10201e6)

        &&

        (((stable_truth_leptons[0].mc_pt > 20e3) && (stable_truth_leptons[1].mc_pt > 27e3))
        ||
        ((stable_truth_leptons[1].mc_pt > 20e3) && (stable_truth_leptons[0].mc_pt > 27e3)))
    );
}

RVec<float> passed_truth_leptons_pt_func(const RVec<TruthParticle>& passed_truth_leptons)
{
    RVec<float> pt;
    pt.reserve(passed_truth_leptons.size());
    for (auto &i: passed_truth_leptons)
    {
        pt.push_back(i.mc_pt/1e3);
    }
    return pt;
}

RVec<float> passed_truth_leptons_eta_func(RVec<TruthParticle>& passed_truth_leptons)
{
    RVec<float> eta;
    eta.reserve(passed_truth_leptons.size());
    for (auto &i: passed_truth_leptons)
    {
        eta.push_back(i.mc_eta);
    }
    return eta;
}

RVec<Photon> selected_photons_func(RVec<Photon> photons)
{
    photons.erase(std::remove_if(photons.begin(),photons.end(),
    [](Photon& x)
    {
        return ((!x.photon_id)|| (x.photon_pt<=10000) || (x.photon_eta >= 2.37)
                || (abs(x.photon_eta)>1.37 && abs(x.photon_eta)<1.52));
        
    }), photons.end());
    
    return photons;
}

RVec<float> selected_photons_pt_func(const RVec<Photon>& selected_photons)
{
    RVec<float> pt;
    pt.reserve(selected_photons.size());
    for (auto &i: selected_photons)
    {
        pt.push_back(i.photon_pt/1e3);
    }
    return pt;
}

RVec<float> selected_photons_eta_func(const RVec<Photon>& selected_photons)
{
    RVec<float> eta;
    eta.reserve(selected_photons.size());
    for (auto &i: selected_photons)
    {
        eta.push_back(i.photon_eta);
    }
    return eta;
}

RVec<float> passed_selected_photons_pt_func(const RVec<Photon>& passed_selected_photons)
{
    RVec<float> pt;
    pt.reserve(passed_selected_photons.size());
    for (auto &i: passed_selected_photons)
    {
        pt.push_back(i.photon_pt/1e3);
    }
    return pt;
}

RVec<float> passed_selected_photons_eta_func(const RVec<Photon>& passed_selected_photons)
{
    RVec<float> eta;
    eta.reserve(passed_selected_photons.size());
    for (auto &i: passed_selected_photons)
    {
        eta.push_back(i.photon_eta);
    }
    return eta;
}

""")


if __name__ == "__main__":
    start = time()
    files = ("/Users/edwardfinkelstein/ATLAS_axion/user.kschmied.28655874._000025.LGNTuple.root",)

    df = MakeRDF(files, 8)
    
    stable_truth_photons = df.Define("stable_truth_photons","stable_truth_photons_func(truth_particles)")
    
    stable_truth_photons_pt = stable_truth_photons.Define("stable_truth_photons_pt", "stable_truth_photons_pt_func(stable_truth_photons)")
    
    stable_truth_photons_eta = stable_truth_photons.Define("stable_truth_photons_eta","stable_truth_photons_eta_func(stable_truth_photons)")
    
    stable_truth_leptons = df.Define("stable_truth_leptons","stable_truth_leptons_func(truth_particles)")
    
    passed_truth_leptons = stable_truth_leptons.Filter("passed_truth_leptons_func(stable_truth_leptons)")
    
    passed_truth_leptons_pt = passed_truth_leptons.Define("passed_truth_leptons_pt", "passed_truth_leptons_pt_func(stable_truth_leptons)")
    
    passed_truth_leptons_eta = passed_truth_leptons.Define("passed_truth_leptons_eta","passed_truth_leptons_eta_func(stable_truth_leptons)")
    
    selected_photons = df.Define("selected_photons","selected_photons_func(photons)")
    
    selected_photons_pt = selected_photons.Define("selected_photons_pt","selected_photons_pt_func(selected_photons)")

    selected_photons_eta = selected_photons.Define("selected_photons_eta","selected_photons_eta_func(selected_photons)")
    
    passed_selected_photons = selected_photons.Filter("selected_photons.size()==2")
    
    passed_selected_photons_pt = passed_selected_photons.Define("passed_selected_photons_pt","passed_selected_photons_pt_func(selected_photons)")
    
    passed_selected_photons_eta = passed_selected_photons.Define("passed_selected_photons_eta","passed_selected_photons_eta_func(selected_photons)")
    
    nEntriesAfterCuts = passed_selected_photons_eta.Count()
    print(f"# events for nominal = {nEntriesAfterCuts.GetValue()}")
    
    histos = (
        stable_truth_photons_pt.Histo1D(("stable_truth_photons_pt", "histTitle", 20, 0, 200), "stable_truth_photons_pt"),
        stable_truth_photons_eta.Histo1D(("stable_truth_photons_eta", "histTitle", 40, -6, 6), "stable_truth_photons_eta"),
        stable_truth_photons_pt.Histo1D(("stable_truth_photons_pt_tight", "histTitle", 100, 0, 100), "stable_truth_photons_pt"),
        passed_truth_leptons_pt.Histo1D(("passed_truth_leptons_pt", "histTitle", 20, 0, 200), "passed_truth_leptons_pt"),
        passed_truth_leptons_eta.Histo1D(("passed_truth_leptons_eta", "histTitle", 40, -6, 6), "passed_truth_leptons_eta"),
        passed_truth_leptons_pt.Histo1D(("passed_truth_leptons_pt_tight", "histTitle", 100, 0, 100), "passed_truth_leptons_pt"),
        selected_photons_pt.Histo1D(("selected_photons_pt", "histTitle", 20, 0, 200), "selected_photons_pt"),
        selected_photons_pt.Histo1D(("selected_photons_pt_tight", "histTitle", 100, 0, 100), "selected_photons_pt"),
        selected_photons_eta.Histo1D(("selected_photons_eta", "histTitle", 40, -6, 6), "selected_photons_eta"),
        passed_selected_photons_pt.Histo1D(("passed_selected_photons_pt_tight", "histTitle", 100, 0, 100), "passed_selected_photons_pt"),
        passed_selected_photons_pt.Histo1D(("passed_selected_photons_pt", "histTitle", 20, 0, 200), "passed_selected_photons_pt"),
        passed_selected_photons_eta.Histo1D(("passed_selected_photons_eta", "histTitle", 40, -6, 6), "passed_selected_photons_eta"),
    )
    
    for h in histos:
        c1 = ROOT.TCanvas("","",800, 700)
        h.Draw("same")
        c1.SaveAs(h.GetName()+".png");
    
    end = time()
    print(f"Time taken = {end-start:0.3f} seconds")
