#Mandatory: List of processes
processList = {
    'p8_ee_ZZ_ecm240':{'chunks':4}, 
#    'p8_ee_WW_ecm240':{'chunks':4}, 
#   'p8_ee_ZH_ecm240':{'chunks':10}, 
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "/eos/user/c/cepeda/FCCeeTests/Taus/test1"

#Optional: analysisName, default is ""
#analysisName = "My Analysis"

#Optional running on HTCondor, default is False
#runBatch    = False

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#Optional test file
testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

#Optional
nCPUS       = 8
runBatch    = False
#batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"



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

               #build pseudo jets with the RP, using the interface that takes px,py,pz,m for better
               #handling of rounding errors
               .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_m)")
               #.Define("pseudo_jets2",    "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

               #run jet clustering with all reconstructed particles. kt_algorithm, R=0.5, exclusive clustering, exactly 4 jets, E0-scheme
               .Define("FCCAnalysesJets_kt", "JetClustering::clustering_kt(0.5, 2, 4, 0, 10)(pseudo_jets)")
               #get the jets out of the struct
               .Define("jets_kt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_kt)")
               #get the jets constituents out of the struct
               .Define("jetconstituents_kt","JetClusteringUtils::get_constituents(FCCAnalysesJets_kt)")
               #get some variables
               .Define("jets_kt_e",        "JetClusteringUtils::get_e(jets_kt)")
               .Define("jets_kt_px",        "JetClusteringUtils::get_px(jets_kt)")
               .Define("jets_kt_py",        "JetClusteringUtils::get_py(jets_kt)")
               .Define("jets_kt_pz",        "JetClusteringUtils::get_pz(jets_kt)")
               .Define("jets_kt_m",        "JetClusteringUtils::get_m(jets_kt)")

               #run jet clustering with all reconstructed particles. ee_genkt_algorithm, R=0.5, inclusive clustering, E-scheme
               .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(0.5, 0, 0, 0, 0, -1)(pseudo_jets)")
               #get the jets out of the struct
               .Define("jets_ee_genkt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")
               #get the jets constituents out of the struct
               .Define("jetconstituents_ee_genkt","JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")
               #get some variables
               .Define("jets_ee_genkt_px",        "JetClusteringUtils::get_px(jets_ee_genkt)")
               .Define("jets_ee_genkt_py",        "JetClusteringUtils::get_py(jets_ee_genkt)")
               .Define("jets_ee_genkt_pz",        "JetClusteringUtils::get_pz(jets_ee_genkt)")
               .Define("jets_ee_genkt_m",        "JetClusteringUtils::get_m(jets_ee_genkt)")


#               #run jet clustering with all reconstructed particles. valencia_algorithm, R=0.5, inclusive clustering, E-scheme
#               .Define("FCCAnalysesJets_valencia", "JetClustering::clustering_valencia(0.5, 1, 6, 0, 0, 1., 1.)(pseudo_jets)")
#
#               #get the jets out of the struct
#               .Define("jets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_valencia)")
#               #get the jets constituents out of the struct
#               .Define("jetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesJets_valencia)")
#               #get some variables
#               .Define("jets_valencia_px",        "JetClusteringUtils::get_px(jets_valencia)")
#               .Define("jets_valencia_py",        "JetClusteringUtils::get_py(jets_valencia)")
#               .Define("jets_valencia_pz",        "JetClusteringUtils::get_pz(jets_valencia)")
#               .Define("jets_valencia_m",        "JetClusteringUtils::get_m(jets_valencia)")
#
#               #run jet clustering with all reconstructed particles. jade_algorithm, R=0.5, exclusive clustering, exactly 4 jets, sorted by E, E0-scheme
#               .Define("FCCAnalysesJets_jade", "JetClustering::clustering_jade(0.5, 2, 4, 1, 10)(pseudo_jets)")
#
#               #get the jets out of the struct
#               .Define("jets_jade",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_jade)")
#               #get the jets constituents out of the struct
#               .Define("jetconstituents_jade","JetClusteringUtils::get_constituents(FCCAnalysesJets_jade)")
#               #get some variables
#               .Define("jets_jade_px",        "JetClusteringUtils::get_px(jets_jade)")
#               .Define("jets_jade_py",        "JetClusteringUtils::get_py(jets_jade)")
#               .Define("jets_jade_pz",        "JetClusteringUtils::get_pz(jets_jade)")
#               .Define("jets_jade_flavour",   "JetTaggingUtils::get_flavour(jets_jade, Particle)")
#               .Define("jets_jade_btag",      "JetTaggingUtils::get_btag(jets_jade_flavour, 0.80)")
#               .Define("jets_jade_btag_true", "JetTaggingUtils::get_btag(jets_jade_flavour, 1.0)")
#               .Define("jets_jade_ctag",      "JetTaggingUtils::get_ctag(jets_jade_flavour, 0.10)")
#               .Define("jets_jade_ctag_true",      "JetTaggingUtils::get_ctag(jets_jade_flavour, 1.0)")
#               .Define("jets_jade_m",        "JetClusteringUtils::get_m(jets_jade)")
#
#               .Define("JET_btag",       "ReconstructedParticle::getJet_btag(Jet3, ParticleIDs, ParticleIDs_0)")
#               .Define("EVT_nbtag",      "ReconstructedParticle::getJet_ntags(JET_btag)")
#
#               .Define('EVT_thrust',     'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
#               .Define('RP_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
#               .Define('EVT_thrust_x',   "EVT_thrust.at(0)")
#               .Define('EVT_thrust_y',   "EVT_thrust.at(1)")
#               .Define('EVT_thrust_z',   "EVT_thrust.at(2)")
#               .Define('EVT_thrust_val', "EVT_thrust.at(3)")
#
#               .Define('EVT_sphericity',     'Algorithms::minimize_sphericity("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
#               .Define('EVT_sphericity_x',   "EVT_sphericity.at(0)")
#               .Define('EVT_sphericity_y',   "EVT_sphericity.at(1)")
#               .Define('EVT_sphericity_z',   "EVT_sphericity.at(2)")
#               .Define('EVT_sphericity_val', "EVT_sphericity.at(3)")
#               .Define('RP_sphericityangle', 'Algorithms::getAxisCosTheta(EVT_sphericity, RP_px, RP_py, RP_pz)')
#
#               .Define('RP_hemis0_mass',   "Algorithms::getAxisMass(0)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
#               .Define('RP_hemis1_mass',   "Algorithms::getAxisMass(1)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
#
#               .Define("RP_total_mass",    "Algorithms::getMass(ReconstructedParticles)")
#
               # Get Electrons
               .Alias("Electron_indices","Electron#0.index")
               .Define("Electron_px","Take(ReconstructedParticles.momentum.x,Electron_indices)")
               .Define("Electron_py","Take(ReconstructedParticles.momentum.y,Electron_indices)")
               .Define("Electron_pz","Take(ReconstructedParticles.momentum.z,Electron_indices)")
               .Define("Electron_mass","Take(ReconstructedParticles.mass,Electron_indices)")
               .Define("Electron_charge","Take(ReconstructedParticles.charge,Electron_indices)")
               .Define("NElectron","(UInt_t)Electron_charge.size()")
               .Define("Electron_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Electron_px, Electron_py, Electron_pz, Electron_mass)")

               # Get Muons
               .Alias("Muon_indices","Muon#0.index")
               .Define("Muon_px","Take(ReconstructedParticles.momentum.x,Muon_indices)")
               .Define("Muon_py","Take(ReconstructedParticles.momentum.y,Muon_indices)")
               .Define("Muon_pz","Take(ReconstructedParticles.momentum.z,Muon_indices)")
               .Define("Muon_mass","Take(ReconstructedParticles.mass,Muon_indices)")
               .Define("Muon_charge","Take(ReconstructedParticles.charge,Muon_indices)")
               .Define("NMuon","(UInt_t)Muon_charge.size()")
               .Define("Muon_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Muon_px, Muon_py, Muon_pz, Muon_mass)")

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
               .Define("Photon_px","Take(ReconstructedParticles.momentum.x,Photon_indices)")
               .Define("Photon_py","Take(ReconstructedParticles.momentum.y,Photon_indices)")
               .Define("Photon_pz","Take(ReconstructedParticles.momentum.z,Photon_indices)")
               .Define("NPhoton","(UInt_t)Photon_px.size()")
               .Define("Photon_mass","Take(ReconstructedParticles.mass,Photon_indices)")
               .Define("Photon_p4","ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > >(Photon_px, Photon_py, Photon_pz, Photon_mass)") # 0.000

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
#               .Define("GenTau_VisPx", "MCParticle::get_px(GenTau_VisP4)")
#               .Define("GenTau_VisTheta", "MCParticle::get_theta(GenTau_VisP4)")
#               .Define("GenTau_VisPhi", "MCParticle::get_phi(GenTau_VisP4)")
#               .Define("GenTau_VisMass", "MCParticle::get_mass(GenTau_VisP4)")
               
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
#               .Define("Tau_pt","ReconstructedParticle::get_pt(tauReco)")
#               .Define("Tau_theta","ReconstructedParticle::get_theta(tauReco)")
#               .Define("Tau_phi","ReconstructedParticle::get_phi(tauReco)")
               .Define("Tau_px","ReconstructedParticle::get_px(tauReco)")
               .Define("Tau_py","ReconstructedParticle::get_py(tauReco)")
               .Define("Tau_pz","ReconstructedParticle::get_pz(tauReco)")
               
               
               
               
        )
        return df2




    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
#                "RP_px", "RP_py", "RP_pz", "RP_e", "RP_m", "RP_q",
#
#                "JET_btag",
#                "EVT_nbtag",
#
#                "EVT_thrust_x", "EVT_thrust_y", "EVT_thrust_z", "EVT_thrust_val",
#
#                "EVT_sphericity_x", "EVT_sphericity_y", "EVT_sphericity_z", "EVT_sphericity_val",
#
#                "RP_thrustangle",
#                "RP_sphericityangle",
#
#                "RP_hemis0_mass",
#                "RP_hemis1_mass",
#                "RP_total_mass",
#
                "jets_kt_e",
                "jets_kt_px",
                "jets_kt_py",
                "jets_kt_pz",
                "jets_kt_m",
                "jetconstituents_kt",

                "jets_ee_genkt_px",
                "jets_ee_genkt_py",
                "jets_ee_genkt_pz",
                "jets_ee_genkt_m",
                "jetconstituents_ee_genkt",

#                "jets_valencia_px",
#                "jets_valencia_py",
#                "jets_valencia_pz",
#                "jets_valencia_m",
#                "jetconstituents_valencia",
#
#                "jets_jade_px",
#                "jets_jade_py",
#                "jets_jade_pz",
#                "jets_jade_m",
#                "jets_jade_ctag",
#                "jets_jade_ctag_true",
#                "jets_jade_btag",
#                "jets_jade_btag_true",
#                "jetconstituents_jade",
#
                "NMuon","Muon_px","Muon_py","Muon_pz","Muon_mass","Muon_charge",
        #        "Muon_energy","Muon_phi","Muon_theta","Muon_pt",
                "NElectron","Electron_px","Electron_py","Electron_pz","Electron_mass","Electron_charge",
        #        "Electron_energy","Electron_phi","Electron_theta","Electron_pt",
                "NPhoton","Photon_px","Photon_py","Photon_pz",
        #        "Photon_energy","Photon_phi","Photon_theta","Photon_pt",
                "Missing_px","Missing_py","Missing_pz", "Missing_pt",
                "Simulated_sqrts",
                "GenTau_px","GenTau_py","GenTau_pz","GenTau_mass","GenTau_PDG",#"GenTau_parent"# ,"GenTau_pt"

#                "GenTau_VisPx",#  "GenTau_VisPhi", "GenTau_VisTheta" , "GenTau_VisMass"
                  "GenTau_VisP4",

                   "pi0RecoTest_mass",

                 "Tau_px" , "Tau_py" , "Tau_pz" , #"Tau_phi", "Tau_theta", 
                 "Tau_mass", "Tau_type",  "Tau_iso" , "Tau_charge" , "NTau" ,



                ]
        return branchList
