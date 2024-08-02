#Mandatory: List of processes
processList = {
    'p8_ee_ZZ_ecm240':{'chunks':598},  # one per job, since two per job failed?
#    'p8_ee_WW_ecm240':{'chunks':20}, 
#   'p8_ee_ZH_ecm240':{'chunks':20}, 
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter23/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "/eos/user/c/cepeda/FCCeeTests/Taus/fileMichele"

#Optional: analysisName, default is ""
#analysisName = "My Analysis"

#Optional running on HTCondor, default is False
runBatch    = False

#Optional batch queue name when running on HTCondor, default is workday
batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
compGroup = "group_u_CMS.u_zh.users" #"group_u_FCC.local_gen"

#Optional test file
testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"
#testFile  ="root://eospublic.cern.ch//eos/user/c/cepeda/FCCeeTests/p8_ee_ZH_vvtautau.root" 

#Optional
nCPUS       = 8



#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               .Alias("Jet3","Jet#3.index")
               #define the RP px, py, pz and e
               .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")
               .Define("RP_m",           "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_q",           "ReconstructedParticle::get_charge(ReconstructedParticles)")

               #select muons on pT
#               .Define("selected_muons", "ReconstructedParticle::sel_pt(5.)(muons)")
               # create branch with muon transverse momentum
#               .Define("zed_leptonic",         "ReconstructedParticle::resonanceBuilder(91)(selected_muons)")
               # create branch with zed transverse momenta
#               .Define("zed_leptonic_recoil",  "ReconstructedParticle::recoilBuilder(240)(zed_leptonic)")


               # Get Electrons
               .Alias("Electron_indices","Electron#0.index")
               .Define("ElectronColl",  "ReconstructedParticle::get(Electron_indices, ReconstructedParticles)")
               .Define("Electron_px","ReconstructedParticle::get_px(ElectronColl)")
               .Define("Electron_py","ReconstructedParticle::get_py(ElectronColl)")
               .Define("Electron_pz","ReconstructedParticle::get_pz(ElectronColl)")
               .Define("Electron_pt","ReconstructedParticle::get_pt(ElectronColl)")
               .Define("Electron_energy","ReconstructedParticle::get_e(ElectronColl)")
               .Define("Electron_theta","ReconstructedParticle::get_theta(ElectronColl)")
               .Define("Electron_phi","ReconstructedParticle::get_phi(ElectronColl)")
               .Define("Electron_mass","ReconstructedParticle::get_mass(ElectronColl)")
               .Define("Electron_y","ReconstructedParticle::get_y(ElectronColl)")
               .Define("Electron_p","ReconstructedParticle::get_p(ElectronColl)")
               .Define("Electron_charge","ReconstructedParticle::get_charge(ElectronColl)")
               .Define("NElectron","(UInt_t)Electron_charge.size()")
               .Define("Electron_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Electron_px, Electron_py, Electron_pz, Electron_mass)")

               # Get Muons
               .Alias("Muon_indices","Muon#0.index")
               .Define("MuonColl",  "ReconstructedParticle::get(Muon_indices, ReconstructedParticles)")
               .Define("Muon_px","ReconstructedParticle::get_px(MuonColl)")
               .Define("Muon_py","ReconstructedParticle::get_py(MuonColl)")
               .Define("Muon_pz","ReconstructedParticle::get_pz(MuonColl)")
               .Define("Muon_pt","ReconstructedParticle::get_pt(MuonColl)")
               .Define("Muon_energy","ReconstructedParticle::get_e(MuonColl)")
               .Define("Muon_theta","ReconstructedParticle::get_theta(MuonColl)")
               .Define("Muon_phi","ReconstructedParticle::get_phi(MuonColl)")
               .Define("Muon_mass","ReconstructedParticle::get_mass(MuonColl)")
               .Define("Muon_y","ReconstructedParticle::get_y(MuonColl)")
               .Define("Muon_p","ReconstructedParticle::get_p(MuonColl)")
               .Define("Muon_charge","ReconstructedParticle::get_charge(MuonColl)")
               .Define("NMuon","(UInt_t)Muon_charge.size()")
               .Define("Muon_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Muon_px, Muon_py, Muon_pz, Muon_mass)")

               # Get Jets
               .Alias("Jet_indices","Muon#0.index")
               .Define("JetColl",  "ReconstructedParticle::get(Jet_indices, ReconstructedParticles)")
               .Define("Jet_px","ReconstructedParticle::get_px(JetColl)")
               .Define("Jet_py","ReconstructedParticle::get_py(JetColl)")
               .Define("Jet_pz","ReconstructedParticle::get_pz(JetColl)")
               .Define("Jet_pt","ReconstructedParticle::get_pt(JetColl)")
               .Define("Jet_energy","ReconstructedParticle::get_e(JetColl)")
               .Define("Jet_theta","ReconstructedParticle::get_theta(JetColl)")
               .Define("Jet_phi","ReconstructedParticle::get_phi(JetColl)")
               .Define("Jet_mass","ReconstructedParticle::get_mass(JetColl)")
               .Define("Jet_y","ReconstructedParticle::get_y(JetColl)")
               .Define("Jet_p","ReconstructedParticle::get_p(JetColl)")
               .Define("Jet_charge","ReconstructedParticle::get_charge(JetColl)")
               .Define("NJet","(UInt_t)Jet_charge.size()")
               .Define("Jet_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Jet_px, Jet_py, Jet_pz, Jet_mass)")

               # Get Missing momentum
               .Define("Missing_px","MissingET.momentum.x[0]")
               .Define("Missing_py","MissingET.momentum.y[0]")
               .Define("Missing_pz","MissingET.momentum.z[0]")
               .Define("Missing_pt","sqrt(Missing_px*Missing_px+Missing_py*Missing_py)")

               # Center-of-mass energy simulated
               .Define("p4_eplus","ROOT::Math::PxPyPzMVector(Particle.momentum.x[0],Particle.momentum.y[0],Particle.momentum.z[0],Particle.mass[0])")
               .Define("p4_eminus","ROOT::Math::PxPyPzMVector(Particle.momentum.x[1],Particle.momentum.y[1],Particle.momentum.z[1],Particle.mass[1])")
               .Define("Simulated_sqrts","(p4_eplus+p4_eminus).M()")
               
               # Get Photons
               .Alias("Photon_indices","Photon#0.index")
               .Define("PhotonColl",  "ReconstructedParticle::get(Photon_indices, ReconstructedParticles)")
               .Define("Photon_px","ReconstructedParticle::get_px(PhotonColl)")
               .Define("Photon_py","ReconstructedParticle::get_py(PhotonColl)")
               .Define("Photon_pz","ReconstructedParticle::get_pz(PhotonColl)")
               .Define("Photon_pt","ReconstructedParticle::get_pt(PhotonColl)")
               .Define("Photon_energy","ReconstructedParticle::get_e(PhotonColl)")
               .Define("Photon_theta","ReconstructedParticle::get_theta(PhotonColl)")
               .Define("Photon_phi","ReconstructedParticle::get_phi(PhotonColl)")
               .Define("Photon_mass","ReconstructedParticle::get_mass(PhotonColl)")
               .Define("Photon_y","ReconstructedParticle::get_y(PhotonColl)")
               .Define("Photon_p","ReconstructedParticle::get_p(PhotonColl)")
               .Define("Photon_charge","ReconstructedParticle::get_charge(PhotonColl)")
               .Define("NPhoton","(UInt_t)Photon_charge.size()")
               .Define("Photon_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Photon_px, Photon_py, Photon_pz, Photon_mass)")

               # Generator information (hard-scattering products)
               .Define("Particle_indices","Nonzero(Particle.PDG)")
               .Alias("Parent_indices","Particle#0.index")
               .Alias("Daughter_indices","Particle#1.index")
               .Define("Particle_parent","Take(Parent_indices,Particle.parents_begin)")
               .Define("Particle_daugther","Take(Parent_indices,Particle.daughters_begin)")
               .Define("Particle_parent_PDG","Take(Particle.PDG,Particle_parent)")
               .Define("Particle_indices_Select","Particle_indices[  abs(Particle.PDG)==15 ]")
               
               .Define("GenTau_PDG","Take(Particle.PDG,Particle_indices_Select)")
               .Define("GenTau_px","Take(Particle.momentum.x,Particle_indices_Select)")
               .Define("GenTau_py","Take(Particle.momentum.y,Particle_indices_Select)")
               .Define("GenTau_pz","Take(Particle.momentum.z,Particle_indices_Select)")
               .Define("GenTau_mass","Take(Particle.mass,Particle_indices_Select)")

               .Define("GenTau_VisP4", "myUtils::VisGenP4(Particle_indices_Select, Particle, Daughter_indices,0)")
               #.Define("GenTau_VisPx", "GenTau_VisP4.X")
               #.Define("GenTau_VisPy", "GenTau_VisP4.Y")
               #.Define("GenTau_VisPz", "GenTau_VisP4.Z")
               #.Define("GenTau_VisPt", "GenTau_VisP4.T")
               #.Define("GenTau_VisTheta", "GenTau_VisP4.Theta")
               #.Define("GenTau_VisPhi", "GenTau_VisP4.Phi")
               #.Define("GenTau_VisMass", "GenTau_VisP4.Mass")
               .Define("NGenTau","(UInt_t)GenTau_VisP4.size()")
 
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

               .Define("pi0RecoTest", "myUtils::MakePi0Test(ReconstructedParticles,0.5)")
               .Define("pi0RecoTest_mass","ReconstructedParticle::get_mass(pi0RecoTest)")
               
               .Define("tauReco", "myUtils::TauRecoTest(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle,pi0RecoTest)")
#               .Define("SortedTaus","SortByIso(tauRecoTest)") g It would be better to sort the RecoParticles, but then I need to sort the indices
               .Define("NTau","(UInt_t)tauReco.size()")
               .Define("Tau_mass","ReconstructedParticle::get_mass(tauReco)")
#               .Define("Tau_p4","ReconstructedParticle::get_tlv(tauReco)")
               .Define("Tau_type","ReconstructedParticle::get_type(tauReco)")
               .Define("Tau_iso","ReconstructedParticle::get_PID(tauReco)")
               .Define("Tau_charge","ReconstructedParticle::get_charge(tauReco)")
               .Define("Tau_pt","ReconstructedParticle::get_pt(tauReco)")
               .Define("Tau_theta","ReconstructedParticle::get_theta(tauReco)")
               .Define("Tau_phi","ReconstructedParticle::get_phi(tauReco)")
               .Define("Tau_px","ReconstructedParticle::get_px(tauReco)")
               .Define("Tau_py","ReconstructedParticle::get_py(tauReco)")
               .Define("Tau_pz","ReconstructedParticle::get_pz(tauReco)")
               .Define("Tau_p","ReconstructedParticle::get_p(tauReco)")
               .Define("Tau_y","ReconstructedParticle::get_y(tauReco)")
               .Define("Tau_energy","ReconstructedParticle::get_e(tauReco)")
               
               
        )
        return df2




    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                "NElectron","Electron_px","Electron_py","Electron_pz","Electron_mass","Electron_charge",
                "Electron_energy","Electron_phi","Electron_theta","Electron_pt",
                "Electron_y","Electron_p",  

                "NMuon","Muon_px","Muon_py","Muon_pz","Muon_mass","Muon_charge",
                "Muon_energy","Muon_phi","Muon_theta","Muon_pt",
                "Muon_y","Muon_p",

                "NJet","Jet_px","Jet_py","Jet_pz","Jet_mass","Jet_charge",
                "Jet_energy","Jet_phi","Jet_theta","Jet_pt",
                "Jet_y","Jet_p",

                "NPhoton","Photon_px","Photon_py","Photon_pz","Photon_mass","Photon_charge",
                "Photon_energy","Photon_phi","Photon_theta","Photon_pt",
                "Photon_y","Photon_p",

                "NTau","Tau_px","Tau_py","Tau_pz","Tau_mass","Tau_charge",
                "Tau_energy","Tau_phi","Tau_theta","Tau_pt",
                "Tau_y","Tau_p",
                "Tau_type",  "Tau_iso" ,

                "Missing_px","Missing_py","Missing_pz", "Missing_pt",
                "Simulated_sqrts",
                "GenTau_px","GenTau_py","GenTau_pz","GenTau_mass","GenTau_PDG",#"GenTau_parent"# ,
                #"GenTau_pt",

                #"GenTau_VisPx", "GenTau_VisPy", "GenTau_VisPz", "GenTau_VisPt","GenTau_VisPhi", "GenTau_VisTheta" , "GenTau_VisMass",
                  "GenTau_VisP4","NGenTau",

                   "pi0RecoTest_mass",

                ]
        return branchList
