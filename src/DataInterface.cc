#include "EventDisplay/inc/DataInterface.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include <sstream>
#include <iomanip>

using namespace mu2e;
namespace REX = ROOT::Experimental;

std::string drawfilename("EventDisplay/config/drawutils.txt");
SimpleConfig drawconfig(drawfilename);


/*
 * Adds reconstructed CaloDigi data products to the REve visualization scene.
 * Digis are visualized as points (crystal center) and 3D boxes.
 * Elements are colored based on the t0 time of the digitization pulse.
*/
void DataInterface::AddCaloDigis(REX::REveManager *&eveMng, bool firstLoop_, 
                                 std::tuple<std::vector<std::string>, 
                                 std::vector<const CaloDigiCollection*>> calodigi_tuple, 
                                 REX::REveElement* &scene){
    /*if(!firstLoop_){
        scene->DestroyElements();;
    }*/
    std::cout << "[DataInterface] AddCaloDigi: Restoring Original Z-Position Logic" << std::endl;
    std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(calodigi_tuple);
    std::vector<std::string> names = std::get<0>(calodigi_tuple);

    // Find the Global Time Range (min/max t0)
    double max_t0 = -1e6;
    double min_t0 = 1e6;
    bool found_digis = false;

    for(const auto* calodigicol : calodigi_list){
        if(calodigicol && !calodigicol->empty()){
            found_digis = true;
            for(const auto& digi : *calodigicol){
                double t0 = digi.t0();
                if(t0 > max_t0) max_t0 = t0;
                if(t0 < min_t0) min_t0 = t0;
            }
        }
    }

    if (!found_digis) {
        std::cout << "[DataInterface] No valid CaloDigis found." << std::endl;
        return;
    }
    
    if (max_t0 == min_t0) {
        max_t0 += 1.0; 
    }

    // Define the Color Gradient (Palette)
    const Int_t NRGBs = 5;
    const Int_t NCont = 255; 
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 }; 
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 }; 
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 }; 
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 }; 
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Visualization Loop
    for(unsigned int j = 0; j < calodigi_list.size(); j++){

        const CaloDigiCollection* calodigicol = calodigi_list[j];
        
        if(calodigicol->size() != 0){
            if(!firstLoop_){
                scene->DestroyElements();;
            }
            
            mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
            GeomHandle<DetectorSystem> det;
            
            auto allcryhits = new REX::REveCompound(("CaloDigiCollection_" + names[j]).c_str(), 
                                                   ("CaloDigi Collection: " + names[j]).c_str(), 1);

            for(unsigned int i = 0; i < calodigicol->size(); i++){
                mu2e::CaloDigi const &digi = (*calodigicol)[i];
                int sipmID = digi.SiPMID();
                int cryID = sipmID / 2;

                Crystal const &crystal = cal.crystal(cryID);
                
                // Calculate Color based on t0
                Color_t color;
                double normalized_t0 = (digi.t0() - min_t0) / (max_t0 - min_t0);
                int colorIdx = static_cast<int>(normalized_t0 * (NCont - 1));
                color = gStyle->GetColorPalette(colorIdx); 

                // Geometry and Position (Matching Original Logic)
                double crystalXLen = pointmmTocm(crystal.size().x());
                double crystalYLen = pointmmTocm(crystal.size().y());
                double crystalZLen = pointmmTocm(crystal.size().z());

                double zpos = 0;
                double diskID = 0;
                if(cryID < 674) {
                    zpos = 186.53;
                    diskID = 0;
                }
                if(cryID >= 674) {
                    zpos = 256.53;
                    diskID = 1;
                }
                // Convert fixed zpos to cm for REve
                double fixed_zpos_cm = zpos;//pointmmTocm(zpos); 

                // Get crystal position in its local Mu2e disk frame (still in mm in Hep3Vector)
                CLHEP::Hep3Vector crystalPos_local_mm = cal.geomUtil().mu2eToDisk(diskID, crystal.position());
                
                // Label
                std::string label = " CaloDigi Instance = " + names[j] + '\n'
                                  + " Crystal ID = " + std::to_string(cryID) + '\n'
                                  + " SiPM. = "+std::to_string(sipmID)+ '\n'
                                  + " t0 = "+std::to_string(digi.t0())+" ns " + '\n'
                                  + " peakPos = "+std::to_string(digi.peakpos());
                
                // A. Draw Crystal Center (Point Set)
                auto ps1 = new REX::REvePointSet(label.c_str(), label.c_str(), 0);
                auto ps2 = new REX::REvePointSet(label.c_str(), label.c_str(), 0);
                
                // The point sets use the fixed global Z position.
                if(diskID == 0)
                    ps1->SetNextPoint(pointmmTocm(crystalPos_local_mm.x()), pointmmTocm(crystalPos_local_mm.y()), fixed_zpos_cm );
                if(diskID == 1) 
                    ps2->SetNextPoint(pointmmTocm(crystalPos_local_mm.x()), pointmmTocm(crystalPos_local_mm.y()), fixed_zpos_cm );

                ps1->SetMarkerColor(color);
                ps1->SetMarkerStyle(DataInterface::mstyle);
                ps1->SetMarkerSize(DataInterface::msize);

                ps2->SetMarkerColor(color);
                ps2->SetMarkerStyle(DataInterface::mstyle);
                ps2->SetMarkerSize(DataInterface::msize);

                scene->AddElement(ps1);
                scene->AddElement(ps2);

                // B. Draw Crystal Volume (REveBox)
                std::string crytitle = "Crystal ID = " + std::to_string(cryID) + '\n'
                                     + " Time = " + std::to_string(digi.t0()) + " ns ";
                
                auto b = new REX::REveBox(crytitle.c_str(), crytitle.c_str());
                b->SetMainColor(color);
                
                double width = crystalXLen / 2.0;
                double height = crystalYLen / 2.0;
                double thickness = crystalZLen / 2.0;
                
                // Z-Terms (All in cm):
                double Z_local_cm = pointmmTocm(crystalPos_local_mm.z());
                
                // The full Z offset term from your original logic: 2*thickness + fixed_zpos_cm
                double Z_offset = 2.0 * thickness + fixed_zpos_cm;

                // Front Face
                b->SetVertex(0, pointmmTocm(crystalPos_local_mm.x()) - width, pointmmTocm(crystalPos_local_mm.y()) - height, Z_local_cm - thickness + Z_offset);
                b->SetVertex(1, pointmmTocm(crystalPos_local_mm.x()) - width, pointmmTocm(crystalPos_local_mm.y()) + height, Z_local_cm - thickness + Z_offset);
                b->SetVertex(2, pointmmTocm(crystalPos_local_mm.x()) + width, pointmmTocm(crystalPos_local_mm.y()) + height, Z_local_cm - thickness + Z_offset);
                b->SetVertex(3, pointmmTocm(crystalPos_local_mm.x()) + width, pointmmTocm(crystalPos_local_mm.y()) - height, Z_local_cm - thickness + Z_offset);
                
                // Back Face (Z_local_cm + thickness + Z_offset)
                b->SetVertex(4, pointmmTocm(crystalPos_local_mm.x()) - width, pointmmTocm(crystalPos_local_mm.y()) - height, Z_local_cm + thickness + Z_offset);
                b->SetVertex(5, pointmmTocm(crystalPos_local_mm.x()) - width, pointmmTocm(crystalPos_local_mm.y()) + height, Z_local_cm + thickness + Z_offset);
                b->SetVertex(6, pointmmTocm(crystalPos_local_mm.x()) + width, pointmmTocm(crystalPos_local_mm.y()) + height, Z_local_cm + thickness + Z_offset);
                b->SetVertex(7, pointmmTocm(crystalPos_local_mm.x()) + width, pointmmTocm(crystalPos_local_mm.y()) - height, Z_local_cm + thickness + Z_offset);

                allcryhits->AddElement(b);
            }
            scene->AddElement(allcryhits);
        }
    }
}


/*
 * Adds reconstructed CaloCluster data products to the REve visualization scene.
 * Crystals in the cluster are visualized as points (crystal center) and 3D boxes.
 * Elements are colored based on the t0 time of the digitization pulse.
 * Crystal lego plot uses energy proportional to length
*/
void DataInterface::AddCaloClusters(REX::REveManager *&eveMng, bool firstLoop_, 
                                    std::tuple<std::vector<std::string>, 
                                    std::vector<const CaloClusterCollection*>> calocluster_tuple, 
                                    REX::REveElement* &scene, bool addCrystalDraw){ // t1 and t2 are now ignored

    /*if(!firstLoop_){
        scene->DestroyElements();;
    }*/
    std::cout << "[DataInterface] AddCaloClusters: Coloring by Proximity to Max-Energy Cluster Time" << std::endl;
    std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(calocluster_tuple);
    std::vector<std::string> names = std::get<0>(calocluster_tuple);

    // Time Difference Color Palette
    // This defines how clusters are colored based on their time separation from t_ref.
    // Near: Red (Prompt) -> Mid: Orange/Yellow -> Far: Green/Blue/Violet (Late/Background)

    // Find the Global Reference Time (t_ref)
    double maxE = 1e-6;
    double t_ref = -1.0; 
    const mu2e::CaloCluster* refCluster = nullptr; // Pointer to the cluster that sets t_ref
    
    // First, iterate over ALL clusters in all collections to find the single most energetic one.
    for(const auto* clustercol : calocluster_list){
        if(clustercol && !clustercol->empty()){
            for(const auto& cluster : *clustercol){
                if(cluster.energyDep() > maxE){
                    maxE = cluster.energyDep();
                    t_ref = cluster.time();
                    refCluster = &cluster;
                }
            }
        }
    }

    // Check if a reference time was found (i.e., at least one cluster exists)
    if(t_ref < 0.0){
        std::cout << "[DataInterface] No valid CaloClusters found to set a reference time." << std::endl;
        return;
    }
    //AddCaloClusterLegend(scene, t_ref); TODO - issue with header file
    // Visualization Loop
    for(unsigned int j = 0; j< calocluster_list.size(); j++){

        const CaloClusterCollection* clustercol = calocluster_list[j];
        
        if(clustercol->size() != 0){
            if(!firstLoop_){
                scene->DestroyElements();;
            }
            
            // Get Geometry Handles
            mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
            GeomHandle<DetectorSystem> det;
            
            // maxE is used here to normalize the crystal box height for the optional Lego plot
            int i = 0;
            // Iterate over the clusters in the current collection
            for(const auto& cluster : *clustercol){
                i++;
                Color_t color = kWhite; // Default
                double markerSize = DataInterface::msize;

                // Calculate time difference
                double time_diff = abs(cluster.time() - t_ref); 

                // Set Color based on Time Proximity to t_ref
                if(&cluster == refCluster){
                    // This is the cluster that set t_ref. Mark it clearly.
                    color = kWhite; 
                    markerSize = DataInterface::msize + 1.0; // Make it larger
                }
                // Clusters close in time to the reference (likely prompt)
                else if(time_diff < 100.0) { // e.g., within +/- 15 ns
                    color = kRed;
                    markerSize = DataInterface::msize;
                }
                else if(time_diff < 250.0) { // e.g., 15-50 ns
                    color = kOrange;
                    markerSize = DataInterface::msize;
                }
                // Clusters far in time (likely background or secondary activity)
                else if(time_diff < 500.0) { // e.g., 50-200 ns
                    color = kYellow;
                    markerSize = DataInterface::msize - 0.5;
                }
                else { // > 200 ns
                    color = kGreen + 2;
                    markerSize = DataInterface::msize - 1.0;
                }

                std::string cluster_energy = std::to_string(cluster.energyDep());
                std::string cluster_time = std::to_string(cluster.time());
                std::string cluster_x = std::to_string(cluster.cog3Vector().x());
                std::string cluster_y = std::to_string(cluster.cog3Vector().y());
                
                CLHEP::Hep3Vector COG(cluster.cog3Vector().x(),cluster.cog3Vector().y(), cluster.cog3Vector().z());
                
                CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(cluster.diskID(),COG);
                CLHEP::Hep3Vector pointInMu2e = det->toMu2e(crystalPos);
                
                std::string cluster_z = std::to_string(abs(pointInMu2e.z()));
                
                std::string label = " Instance = " + names[0] + '\n'
                                  + " E = "+cluster_energy+" MeV " + '\n'
                                  + " Time = "+cluster_time+" ns " + '\n'
                                  + " |dt| from Emax = " + std::to_string(time_diff) + " ns " + '\n'
                                  + " Pos = ("+cluster_x+","+cluster_y+","+cluster_z+") mm";
                
                std::string name = "disk" + std::to_string(cluster.diskID()) + label;
                auto ps1 = new REX::REvePointSet(name, "CaloClusters Disk 1: "+label,0);
                auto ps2 = new REX::REvePointSet(name, "CaloClusters Disk 2: "+label,0);


                // Set Positions, Color, and Size
                if(cluster.diskID() == 0)    
                    ps1->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));
                if(cluster.diskID() == 1)
                    ps2->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));
                
                ps1->SetMarkerColor(color);
                ps1->SetMarkerStyle(DataInterface::mstyle);
                ps1->SetMarkerSize(markerSize);

                ps2->SetMarkerColor(color);
                ps2->SetMarkerStyle(DataInterface::mstyle);
                ps2->SetMarkerSize(markerSize);

                scene->AddElement(ps1);
                scene->AddElement(ps2);

                // Optional: Draw Contributing Crystals (Lego Plot)
                if(addCrystalDraw){
                    auto allcryhits = new REX::REveCompound("CrystalHits for Cluster "+std::to_string(cluster.diskID())+"_"+std::to_string(i),"CrystalHits for Cluster "+std::to_string(cluster.diskID())+"_"+std::to_string(i),1);
                    
                    for(unsigned h =0 ; h < cluster.caloHitsPtrVector().size();h++)    {
                        art::Ptr<CaloHit> crystalhit = cluster.caloHitsPtrVector()[h];
                        int cryID = crystalhit->crystalID();

                        Crystal const &crystal = cal.crystal(cryID);
                        double crystalXLen = pointmmTocm(crystal.size().x());
                        double crystalYLen = pointmmTocm(crystal.size().y());
                        double crystalZLen = pointmmTocm(crystal.size().z());

                        CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(cluster.diskID(),crystal.position()) ;

                        std::string crytitle = "disk"+std::to_string(cal.crystal(crystalhit->crystalID()).diskID()) + " Crystal Hit = " + std::to_string(cryID) + '\n'
                                             + " Energy Dep. = "+std::to_string(crystalhit->energyDep())+" MeV " + '\n'
                                             + " Time = "+std::to_string(crystalhit->time())+" ns ";
                        
                        char const *crytitle_c = crytitle.c_str();
                        auto b = new REX::REveBox(crytitle_c,crytitle_c); 

                        b->SetMainColor(color); // Use the time-proximity color

                        double width = crystalXLen/2;
                        double height = crystalYLen/2;
                        
                        // Z proportional to energy (normalized by the global max E found earlier)
                        double thickness = crystalhit->energyDep()/maxE * crystalZLen/2; 
                        
                        // Calculate Z offset
                        double crystalZOffset = pointmmTocm(crystalPos.z()) + abs(pointmmTocm(pointInMu2e.z())) + crystalZLen/2;
                        
                        // Set the 8 vertices of the REveBox object
                        b->SetVertex(0, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height , crystalZOffset - thickness); 
                        b->SetVertex(1, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height, crystalZOffset - thickness); 
                        b->SetVertex(2, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height ,crystalZOffset - thickness); 
                        b->SetVertex(3, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, crystalZOffset - thickness); 
                        b->SetVertex(4, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,crystalZOffset + thickness); 
                        b->SetVertex(5, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height , crystalZOffset + thickness); 
                        b->SetVertex(6, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height , crystalZOffset + thickness); 
                        b->SetVertex(7,pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, crystalZOffset + thickness); 
                        
                        allcryhits->AddElement(b);
                    }
                    scene->AddElement(allcryhits);
                }
            }
        }
    }
}

/*
Enables the visualization of cluster of hits flagged as background by the FlagBkgHits module. This is work in progress. More features coming soon.
*/
void DataInterface::AddBkgClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const BkgClusterCollection*>> bkgcluster_tuple, REX::REveElement* &scene){
  std::cout<<"BkgClusterCollection "<<std::endl;
  std::vector<const BkgClusterCollection*> bkgcluster_list = std::get<1>(bkgcluster_tuple);
  // std::vector<std::string> names = std::get<0>(bkgcluster_tuple);
  std::cout<<"BkgClusterCollection size = "<<bkgcluster_list.size()<<std::endl;
  for(unsigned int j = 0; j < bkgcluster_list.size(); j++){
    const BkgClusterCollection* bccol = bkgcluster_list[j];
    if(bccol->size() !=0 ){
      // Loop over hits
      for(unsigned int i=0; i< bccol->size(); i++){
        mu2e::BkgCluster const  &bkgcluster= (*bccol)[i];
        int colour = (i+3);
        std::cout<<"BkgCluster ="<<bkgcluster.hits().size()<<std::endl;
        CLHEP::Hep3Vector ClusterPos(pointmmTocm(bkgcluster.pos().x()), pointmmTocm(bkgcluster.pos().y()), pointmmTocm(bkgcluster.pos().z()));
        std::string bctitle = "BkgCluster";
        auto ps1 = new REX::REvePointSet(bctitle, bctitle,0);
        ps1->SetNextPoint(ClusterPos.x(), ClusterPos.y() , ClusterPos.z());
        ps1->SetMarkerColor(colour);
        ps1->SetMarkerStyle(DataInterface::mstyle);
        ps1->SetMarkerSize(DataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }
}

/*
 * Adds reconstructed ComboHits data products to the REve visualization scene.
 * Hits are visualized as points with optional error bars
 * Elements are colored based on the t0 time of the digitization pulse.
*/

//FIXME if firstloop never used remove it
void DataInterface::AddComboHits(REX::REveManager *&eveMng, bool firstLoop_, 
                                 std::tuple<std::vector<std::string>, 
                                 std::vector<const ComboHitCollection*>> combohit_tuple, 
                                 REX::REveElement* &scene, 
                                 bool strawdisplay, bool AddErrorBar_) {
    /*if(!firstLoop_){
        scene->DestroyElements();;
    }*/
    std::vector<const ComboHitCollection*> combohit_list = std::get<1>(combohit_tuple);
    std::vector<std::string> names = std::get<0>(combohit_tuple);

    // Find the Global Time Range (min/max time)
    double max_time = -1e6;
    double min_time = 1e6;
    bool found_hits = false;

    for(const auto* chcol : combohit_list){
        if(chcol && !chcol->empty()){
            found_hits = true;
            for(const auto& hit : *chcol){
                double time = hit.time();
                if(time > max_time) max_time = time;
                if(time < min_time) min_time = time;
            }
        }
    }

    if (!found_hits) {
        std::cout << "[DataInterface] No valid ComboHits found." << std::endl;
        return;
    }
    if (max_time == min_time) {
        max_time += 1.0; 
    }

    // Define the Color Gradient (Palette)
    const Int_t NRGBs = 5;
    const Int_t NCont = 255; 
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 }; 
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 }; 
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 }; 
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 }; 
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Visualization Loop
    for(unsigned int j = 0; j < combohit_list.size(); j++){
        const ComboHitCollection* chcol = combohit_list[j];
        
        if(chcol->size() != 0){
            
            // MASTER COMPOUND for the entire collection
            std::string master_name = "ComboHits_" + names[j];
            auto master_compound = new REX::REveCompound(master_name.c_str(), master_name.c_str(), true);
            master_compound->SetRnrSelf(false); // Only render children (hits/errors)
            
            mu2e::GeomHandle<mu2e::Tracker> tracker;
            const auto& allStraws = tracker->getStraws();
            
            // Loop over hits
            for(unsigned int i = 0; i < chcol->size(); i++){
                mu2e::ComboHit const &hit = (*chcol)[i];
                
                // Calculate Time-Based Color
                Color_t hit_color;
                double normalized_time = (hit.time() - min_time) / (max_time - min_time);
                int colorIdx = static_cast<int>(normalized_time * (NCont - 1));
                hit_color = gStyle->GetColorPalette(colorIdx); 

                // A. Display Hit Straws (Separate element, added directly to scene for broad context)
                if(strawdisplay){
                    int sid = hit._sid.asUint16();
                    
                    if(static_cast<std::size_t>(sid) < allStraws.size() && sid < drawconfig.getInt("maxSID")){
                        CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
                        const mu2e::Straw& s = allStraws[sid];
                        const CLHEP::Hep3Vector& p = s.getMidPoint();
                        const CLHEP::Hep3Vector& d = s.getDirection();
                        double x = p.x();
                        double y = p.y();
                        double z = p.z();
                        double l = s.halfLength();
                        
                        double st = sin(d.theta());
                        double ct = cos(d.theta());
                        double sp = sin(d.phi()+TMath::Pi()/2.0);
                        double cp = cos(d.phi()+TMath::Pi()/2.0);

                        double x1=x+l*st*sp;
                        double y1=y-l*st*cp;
                        double z1=z+l*ct;
                        double x2=x-l*st*sp;
                        double y2=y+l*st*cp;
                        double z2=z-l*ct;

                        std::string strawtitle;
                        int idPlane = s.id().getPlane();
                        int colorid = s.id().getPanel() + idPlane + 1;
                        strawtitle =Form("Straw %i Panel %i Plane %i",s.id().getStraw(), s.id().getPanel(), idPlane);
                        
                        sposi.set(x1, y1, z1);
                        sposf.set(x2, y2, z2);
                        
                        if(sposi.x() != 0){
                            auto strawline = new REX::REveLine("StrawHit", strawtitle, 2);
                            strawline->SetPoint(0, pointmmTocm(sposi.x()), pointmmTocm(sposi.y()), pointmmTocm(sposi.z()));
                            strawline->SetNextPoint(pointmmTocm(sposf.x()), pointmmTocm(sposf.y()), pointmmTocm(sposf.z()));
                            strawline->SetLineWidth(1);
                            strawline->SetLineColor(colorid); 
                            if(strawline->GetSize() != 0) scene->AddElement(strawline);
                        }
                    }
                }

                // B. Add Error Bar (REveLine) Optional
                if(AddErrorBar_){
                    auto const& p = hit.pos();
                    auto w = hit.uDir();
                    auto const& s = hit.wireRes(); // Wire resolution (length of error bar half-segment)
                    
                    // Calculate endpoints of the error bar along the perpendicular direction (w = uDir)
                    double x1 = (p.x()+s*w.x());
                    double x2 = (p.x()-s*w.x());
                    double z1 = (p.z()+s*w.z());
                    double z2 = (p.z()-s*w.z());
                    double y1 = (p.y()+s*w.y());
                    double y2 = (p.y()-s*w.y());

                    // Add a detailed error bar label
                    std::string errorbar_label = std::string("ComboHit Error Bar") + '\n' 
                                               + "Hit Time: " + std::to_string(hit.time()) + " ns" + '\n'
                                               + "Wire Res (half-length): " + std::to_string(s) + " mm" + '\n'
                                               + "Error Bar Endpoints (mm):" + '\n'
                                               + " P1 (" + std::to_string(x1) + ", " + std::to_string(y1) + ", " + std::to_string(z1) + ")" + '\n'
                                               + " P2 (" + std::to_string(x2) + ", " + std::to_string(y2) + ", " + std::to_string(z2) + ")";
                    auto error = new REX::REveLine(errorbar_label.c_str(), errorbar_label.c_str(), 2);
                    error->SetPoint(0, pointmmTocm(x1), pointmmTocm(y1), pointmmTocm(z1));
                    error->SetNextPoint(pointmmTocm(x2), pointmmTocm(y2), pointmmTocm(z2));
                    
                    error->SetLineColor(hit_color); 
                    error->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
                    
                    master_compound->AddElement(error);
                }
                
                // C. Draw ComboHit Position (REvePointSet)
                CLHEP::Hep3Vector HitPos(pointmmTocm(hit.pos().x()), pointmmTocm(hit.pos().y()), pointmmTocm(hit.pos().z()));
                std::string chtitle = "ComboHits tag = "
                    + (names[j]) + '\n'
                    + " position : x " + std::to_string(hit.pos().x()) + '\n'
                    + " y " + std::to_string(hit.pos().y()) + '\n'
                    + " z " + std::to_string(hit.pos().z()) + '\n'
                    + " time :" + std::to_string(hit.time()) + '\n'
                    + " energy dep : " + std::to_string(hit.energyDep()) + "MeV";
                    
                auto ps1 = new REX::REvePointSet(chtitle.c_str(), chtitle.c_str(), 0);
                ps1->SetNextPoint(HitPos.x(), HitPos.y() , HitPos.z());
                
                ps1->SetMarkerColor(hit_color); // Time-based color
                ps1->SetMarkerStyle(DataInterface::mstyle);
                ps1->SetMarkerSize(DataInterface::msize);
                
                master_compound->AddElement(ps1);
            } 
            
            // Add the Master Compound to the Scene
            scene->AddElement(master_compound);
        }
    }
}


/*
 * Adds reconstructed CaloRecoPulse data products to the REve visualization scene.
 * RecoPules are visualized as points bars are visualized as 3D boxes
 * Elements are colored based on the height of the digitization pulse.
*/
void DataInterface::AddCrvInfo(REX::REveManager *&eveMng, bool firstLoop_, 
                               std::tuple<std::vector<std::string>, 
                               std::vector<const CrvRecoPulseCollection*>> crvpulse_tuple, 
                               REX::REveElement* &scene, bool extracted, bool addCrvBars){

    /*if(!firstLoop_){
        scene->DestroyElements();;
    }*/
    std::cout << "[DataInterface::AddCrvInfo() ] - Coloring by Pulse Height (White Start, More Variation)" << std::endl;
    std::vector<const CrvRecoPulseCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
    std::vector<std::string> names = std::get<0>(crvpulse_tuple);
    
    // Geometry Handles
    mu2e::GeomHandle<mu2e::CosmicRayShield> CRS;
    mu2e::GeomHandle<mu2e::DetectorSystem> det;
    
    // Find the Global Pulse Height Range (min/max)
    double max_height = 0.0;
    double min_height = 1e6;
    bool found_pulses = false;

    if(crvpulse_list.size() != 0){
        for(const auto* crvRecoPulse : crvpulse_list){
            if(crvRecoPulse && !crvRecoPulse->empty()){
                found_pulses = true;
                for(const auto& crvpulse : *crvRecoPulse){
                    double height = crvpulse.GetPulseHeight();
                    if(height > max_height) max_height = height;
                    if(height < min_height) min_height = height;
                }
            }
        }
    }

    if (!found_pulses) {
        std::cout << "[ REveDataInterface::AddCrvInfo() ] No valid Crv Reco Pulses found." << std::endl;
        return;
    }
    
    // Handle cases where all pulses have the same height
    if (max_height == min_height) {
        max_height += 1.0; 
    }

    // Define the Color Gradient (Palette)
    const Int_t NRGBs = 7; 
    const Int_t NCont = 255; 

    Double_t stops[NRGBs] = { 0.00, 0.15, 0.30, 0.50, 0.70, 0.85, 1.00 }; 
    Double_t red[NRGBs]   = { 1.00, 0.90, 0.60, 0.00, 0.00, 0.50, 0.80 }; 
    Double_t green[NRGBs] = { 1.00, 0.90, 0.80, 0.60, 0.00, 0.00, 0.10 }; 
    Double_t blue[NRGBs]  = { 1.00, 0.90, 0.80, 0.80, 1.00, 0.80, 0.00 }; 
    
    Int_t palette_start_index = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    
    // Visualization Loop
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
        const CrvRecoPulseCollection* crvRecoPulse = crvpulse_list[i];
        
        if(crvRecoPulse->size() != 0){
            
            std::string temp_title = "Crv Reco Pulses: " + names[i];
            auto ps1 = new REX::REvePointSet(temp_title.c_str(), temp_title.c_str(), 0); 
            auto allcrvbars = new REX::REveCompound(("CrvHitBars_" + std::to_string(i)).c_str(), 
                                                   ("Crv Hit Bars: " + names[i]).c_str(), 1);

            std::string final_ps1_title = temp_title;
            
            for(unsigned int j=0; j< crvRecoPulse->size(); j++){

                mu2e::CrvRecoPulse const &crvpulse = (*crvRecoPulse)[j];
                
                // Calculate Pulse Height Based Color
                Color_t hit_color;
                
                // Define a minimum range (e.g., 1.0 ADC count) to prevent color collapse
                const double kMinRange = 1.0; 

                double pulse_height = crvpulse.GetPulseHeight();
                double actual_range = max_height - min_height;

                // Use the larger of the actual range or the minimum required range (kMinRange)
                double range = (actual_range < kMinRange) ? kMinRange : actual_range;

                // Calculate the normalized height (0.0 to 1.0) using the adjusted range
                double normalized_height = (pulse_height - min_height) / range;

                // Clamp the normalized height between 0 and 1
                if (normalized_height < 0) normalized_height = 0;
                if (normalized_height > 1) normalized_height = 1;

                // Map normalized height (0 to 1) to the color index (0 to NCont-1)
                int colorIdx = palette_start_index + static_cast<int>(normalized_height * (NCont - 1));

                hit_color = TColor::GetColor(colorIdx);


                // Geometry Setup
                const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvpulse.GetScintillatorBarIndex();
                const mu2e::CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
                const mu2e::CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
                
                CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition(); 
                CLHEP::Hep3Vector pointInMu2e = det->toDetector(crvCounterPos); 

                // DETAILED PULSE TITLE (Used for REveBox mouseover)
                std::string pulsetitle = " Crv Pulse Collection: " + names[i] + '\n'
                                       + " SiPM Channel ID: " + std::to_string(crvpulse.GetSiPMNumber()) + '\n'
                                       + " Crv Bar ID: " + std::to_string(crvBarIndex.asInt()) + '\n'
                                       + " Pulse Time: " + std::to_string(crvpulse.GetPulseTime()) + " ns" + '\n'
                                       + " Pulse Height: " + std::to_string(crvpulse.GetPulseHeight()) + " ADC (Color)" + '\n'
                                       + " Pedestal: " + std::to_string(crvpulse.GetPedestal()) + " ADC";
                char const *pulsetitle_c = pulsetitle.c_str();

                final_ps1_title = "Crv Pulses: " + names[i] + " (Last Pulse Details: " + '\n' + pulsetitle;
                
                // A. Draw Crv Bar (REveBox)
                if(addCrvBars){
                   
                    double length = pointmmTocm(crvCounter.getHalfLength());
                    double width = pointmmTocm(crvCounter.getHalfWidth());
                    double height = pointmmTocm(crvCounter.getHalfThickness());
                    
                    double X_cm = pointmmTocm(pointInMu2e.x());
                    double Y_cm = pointmmTocm(pointInMu2e.y());
                    double Z_cm = pointmmTocm(pointInMu2e.z());
                    
                    // Handle geometry errors that cause Seg Faults
                    try {
                        // Check for degenerate dimensions (re-added for robustness)
                        const double kMinDimension = 1e-4; // 1 micrometer (in cm)
                        if (length < kMinDimension || width < kMinDimension || height < kMinDimension) {
                             throw std::runtime_error("Degenerate dimensions detected.");
                        }

                        // Only create and configure the REveBox if the dimensions are safe
                        auto b = new REX::REveBox(pulsetitle_c, pulsetitle_c);
                        b->SetMainColor(hit_color); 
                        b->SetMainTransparency(drawconfig.getInt("Crvtrans"));
                        b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
                        if(!extracted){
                          // Crv D, Crv U orientation (Length along X)
                          if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){ 
                              b->SetVertex(0, X_cm - length, Y_cm - width, Z_cm - height);
                              b->SetVertex(1, X_cm - length, Y_cm - width, Z_cm + height);
                              b->SetVertex(2, X_cm - length, Y_cm + width, Z_cm + height);
                              b->SetVertex(3, X_cm - length, Y_cm + width, Z_cm - height);
                              b->SetVertex(4, X_cm + length, Y_cm - width, Z_cm - height);
                              b->SetVertex(5, X_cm + length, Y_cm - width, Z_cm + height);
                              b->SetVertex(6, X_cm + length, Y_cm + width, Z_cm + height);
                              b->SetVertex(7, X_cm + length, Y_cm + width, Z_cm - height);
                          }
                          // Crv TS orientation (Length along Z)
                          else if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){
                              b->SetVertex(0, X_cm - width, Y_cm - height, Z_cm - length);
                              b->SetVertex(1, X_cm - width, Y_cm + height, Z_cm - length);
                              b->SetVertex(2, X_cm + width, Y_cm + height, Z_cm - length);
                              b->SetVertex(3, X_cm + width, Y_cm - height, Z_cm - length);
                              b->SetVertex(4, X_cm - width, Y_cm - height, Z_cm + length);
                              b->SetVertex(5, X_cm - width, Y_cm + height, Z_cm + length);
                              b->SetVertex(6, X_cm + width, Y_cm + height, Z_cm + length);
                              b->SetVertex(7, X_cm + width, Y_cm - height, Z_cm + length);
                          }
                          // Crv T orientation (Length along X)
                          else if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){ 
                              b->SetVertex(0, X_cm - length, Y_cm - height, Z_cm - width);
                              b->SetVertex(1, X_cm - length, Y_cm + height, Z_cm - width);
                              b->SetVertex(2, X_cm - length, Y_cm + height, Z_cm + width);
                              b->SetVertex(3, X_cm - length, Y_cm - height, Z_cm + width);
                              b->SetVertex(4, X_cm + length, Y_cm - height, Z_cm - width);
                              b->SetVertex(5, X_cm + length, Y_cm + height, Z_cm - width);
                              b->SetVertex(6, X_cm + length, Y_cm + height, Z_cm + width);
                              b->SetVertex(7, X_cm + length, Y_cm - height, Z_cm + width);
                          }
                          // Crv R, Crv L orientation (Length along Y) 
                          else if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 1){
                              b->SetVertex(0, X_cm - height, Y_cm - length, Z_cm - width);
                              b->SetVertex(1, X_cm - height, Y_cm + length, Z_cm - width);
                              b->SetVertex(2, X_cm + height, Y_cm + length, Z_cm - width);
                              b->SetVertex(3, X_cm + height, Y_cm - length, Z_cm - width); 
                              
                              b->SetVertex(4, X_cm - height, Y_cm - length, Z_cm + width);
                              b->SetVertex(5, X_cm - height, Y_cm + length, Z_cm + width);
                              b->SetVertex(6, X_cm + height, Y_cm + length, Z_cm + width);
                              b->SetVertex(7, X_cm + height, Y_cm - length, Z_cm + width);
                          }
                        }
                        if(extracted){
                            if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //T1 (Length along Z)
                                b->SetVertex(0, X_cm - width, Y_cm - height, -1*length + Z_cm);
                                b->SetVertex(1, X_cm - width, Y_cm + height, -1*length + Z_cm);
                                b->SetVertex(2, X_cm + width, Y_cm + height, -1*length + Z_cm);
                                b->SetVertex(3, X_cm + width, Y_cm - height, -1*length + Z_cm);
                                b->SetVertex(4, X_cm - width, Y_cm - height, length + Z_cm);
                                b->SetVertex(5, X_cm - width, Y_cm + height, length + Z_cm);
                                b->SetVertex(6, X_cm + width, Y_cm + height, length + Z_cm);
                                b->SetVertex(7, X_cm + width, Y_cm - height, length + Z_cm);
                            } else { //EX, T2 (Length along X)
                                b->SetVertex(0, X_cm - length, Y_cm - height, Z_cm - width);
                                b->SetVertex(1, X_cm - length, Y_cm + height, Z_cm - width);
                                b->SetVertex(2, X_cm + length, Y_cm + height, Z_cm - width);
                                b->SetVertex(3, X_cm + length, Y_cm - height , Z_cm - width);
                                b->SetVertex(4, X_cm - length , Y_cm - height, Z_cm + width);
                                b->SetVertex(5, X_cm - length, Y_cm + height , Z_cm + width);
                                b->SetVertex(6, X_cm + length, Y_cm + height , Z_cm + width);
                                b->SetVertex(7, X_cm + length , Y_cm - height, Z_cm + width);
                            }
                            scene->AddElement(b); 
                        } else {
                             allcrvbars->AddElement(b); 
                        }
                    } catch (const std::exception& e) {
                        // Print a warning but continue the loop
                        std::cerr << "Crv Geometry Error (Bar ID: " << crvBarIndex.asInt() 
                                  << "). Skipping REveBox creation. Reason: " << e.what() 
                                  << " L=" << length << ", W=" << width << ", H=" << height << std::endl;
                    } catch (...) {
                        // Catch any unknown exception type (e.g., if REve throws a non-standard error)
                        std::cerr << "Crv Geometry Critical Error (Bar ID: " << crvBarIndex.asInt() 
                                  << "). Skipping REveBox creation due to unknown geometry failure." << std::endl;
                    }
                }
                
                // B. Add Reco Pulse Marker (Point Set)
                ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()), pointmmTocm(pointInMu2e.z()));
                
            } // End of individual pulse loop
            
            // Post-Loop: Set the final, detailed title on the PointSet
            ps1->SetTitle(final_ps1_title.c_str());
            
            // Add Compounds/Collections to Scene
            if (!extracted && addCrvBars) {
                scene->AddElement(allcrvbars);
            }
            
            // Draw reco pulse collection (markers)
            ps1->SetMarkerColor(i + 3); 
            ps1->SetMarkerStyle(DataInterface::mstyle);
            ps1->SetMarkerSize(DataInterface::msize);
            
            if(ps1->GetSize() != 0) scene->AddElement(ps1);
        }
    }
}

/*
 * Adds reconstructed CrvCluster data products to the REve visualization scene.
 * Clusters are visualized as points, RecoPulses used to access bar info, bars visualized as boxes in 3D
 * Elements are colored based on the pulse height of the digitization pulse.
*/
void DataInterface::AddCrvClusters(REX::REveManager *&eveMng, bool firstLoop_, 
                                 std::tuple<std::vector<std::string>, 
                                 std::vector<const CrvCoincidenceClusterCollection*>> crvpulse_tuple, 
                                 REX::REveElement* &scene, bool extracted, bool addCrvBars)
{

    std::vector<const CrvCoincidenceClusterCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
    std::vector<std::string> names = std::get<0>(crvpulse_tuple);
    mu2e::GeomHandle<mu2e::CosmicRayShield> CRS;
    mu2e::GeomHandle<mu2e::DetectorSystem> det;
    
    if (crvpulse_list.size() == 0) return;

    // Find the Global Pulse Height Range (min/max) for ALL pulses in ALL clusters
    double max_height = 0.0;
    double min_height = 1e6;
    bool found_pulses = false;

    for (const auto* crvClusters : crvpulse_list) {
        if (crvClusters) {
            for (const auto& crvclu : *crvClusters) {
                for (const auto& crvpulsePtr : crvclu.GetCrvRecoPulses()) {
                    if (crvpulsePtr) {
                        found_pulses = true;
                        double height = crvpulsePtr->GetPulseHeight();
                        if(height > max_height) max_height = height;
                        if(height < min_height) min_height = height;
                    }
                }
            }
        }
    }

    if (!found_pulses) {
        std::cout << "[DataInterface::AddCrvClusters() ] No valid Crv Reco Pulses found in clusters." << std::endl;
        return;
    }

    // Handle cases where all pulses have the same height
    if (max_height == min_height) {
        max_height += 1.0; 
    }
    
    // Define the Color Gradient (Palette)
    const Int_t NRGBs = 7;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.15, 0.30, 0.50, 0.70, 0.85, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.90, 0.60, 0.00, 0.00, 0.50, 0.80 };
    Double_t green[NRGBs] = { 1.00, 0.90, 0.80, 0.60, 0.00, 0.00, 0.10 };
    Double_t blue[NRGBs]  = { 1.00, 0.90, 0.80, 0.80, 1.00, 0.80, 0.00 };
    Int_t palette_start_index = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    const double kMinRange = 1.0; // Define a minimum range to prevent color collapse

    // Visualization Loop
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
        const CrvCoincidenceClusterCollection* crvClusters = crvpulse_list[i];
        
        if (crvClusters->size() == 0) continue;

        std::string bars_title = "Crv Cluster Bars: " + names[i];
        auto allcrvbars_collection = new REX::REveCompound(bars_title.c_str(), bars_title.c_str(), 1); 

        for(unsigned int j=0; j< crvClusters->size(); j++){
            
            mu2e::CrvCoincidenceCluster const &crvclu = (*crvClusters)[j];
            
            // Make title
            std::string crvtitle = "CrvCoincidenceCluster" + std::to_string(j) + " tag : " + names[i] + '\n'
            + " averge hit time = " + std::to_string(crvclu.GetAvgHitTime())+" ns " + '\n'
            + " PEs = " + std::to_string(crvclu.GetPEs());
            auto ps1 = new REX::REvePointSet(crvtitle.c_str(), crvtitle.c_str(), 0);
            
            CLHEP::Hep3Vector pointInMu2e = det->toDetector(crvclu.GetAvgHitPos());
            ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()) , pointmmTocm(pointInMu2e.z()));
            
            std::set<mu2e::CRSScintillatorBarIndex> drawn_bars; 

            for(unsigned h =0 ; h < crvclu.GetCrvRecoPulses().size();h++) {
                
                art::Ptr<mu2e::CrvRecoPulse> crvpulse = crvclu.GetCrvRecoPulses()[h];
                const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvpulse->GetScintillatorBarIndex();

                if (addCrvBars && crvpulse && drawn_bars.count(crvBarIndex) == 0) {
                    
                    drawn_bars.insert(crvBarIndex);
                    
                    const mu2e::CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
                    const mu2e::CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
                    CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
                    
                    // Calculate Pulse Height Based Color
                    double pulse_height = crvpulse->GetPulseHeight();
                    double actual_range = max_height - min_height;
                    double range = (actual_range < kMinRange) ? kMinRange : actual_range;

                    double normalized_height = (pulse_height - min_height) / range;
                    if (normalized_height < 0) normalized_height = 0;
                    if (normalized_height > 1) normalized_height = 1;

                    int colorIdx = palette_start_index + static_cast<int>(normalized_height * (NCont - 1));
                    Color_t hit_color = TColor::GetColor(colorIdx);


                    // Base Geometry Setup
                    CLHEP::Hep3Vector bar_center_mu2e = det->toDetector(crvCounterPos);
                    
                    double X_cm = pointmmTocm(bar_center_mu2e.x());
                    double Y_cm = pointmmTocm(bar_center_mu2e.y());
                    double Z_cm = pointmmTocm(bar_center_mu2e.z());

                    double length = pointmmTocm(crvCounter.getHalfLength());
                    double width = pointmmTocm(crvCounter.getHalfWidth());
                    double height = pointmmTocm(crvCounter.getHalfThickness());

                    std::string pulsetitle = " Crv Bar Hit for tag : "
                                            + names[i] +  '\n'
                                            + "Bar ID: " + std::to_string(crvBarIndex.asInt()) +  '\n'
                                            + "Pulse Height: " + std::to_string(pulse_height) + " ADC (Color)" + '\n'
                                            + "Coincidence start time: " + std::to_string(crvclu.GetStartTime()) +  '\n'
                                            + "Coincidence end time: " + std::to_string(crvclu.GetEndTime());
                    char const *pulsetitle_c = pulsetitle.c_str();

                    auto b = new REX::REveBox(pulsetitle_c, pulsetitle_c);
                    b->SetMainColor(hit_color);
                    b->SetMainTransparency(drawconfig.getInt("Crvtrans"));
                    b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
                    
                    // Geometry Setting
                    if(!extracted){
                        // Crv D, Crv U orientation (Length along X)
                        if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){
                            b->SetVertex(0, X_cm - length, Y_cm - width, Z_cm - height);
                            b->SetVertex(1, X_cm - length, Y_cm - width, Z_cm + height);
                            b->SetVertex(2, X_cm - length, Y_cm + width, Z_cm + height);
                            b->SetVertex(3, X_cm - length, Y_cm + width, Z_cm - height);
                            b->SetVertex(4, X_cm + length, Y_cm - width, Z_cm - height);
                            b->SetVertex(5, X_cm + length, Y_cm - width, Z_cm + height);
                            b->SetVertex(6, X_cm + length, Y_cm + width, Z_cm + height);
                            b->SetVertex(7, X_cm + length, Y_cm + width, Z_cm - height);
                        }
                        // Crv TS orientation (Length along Z)
                        else if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){
                            b->SetVertex(0, X_cm - width, Y_cm - height, Z_cm - length);
                            b->SetVertex(1, X_cm - width, Y_cm + height, Z_cm - length);
                            b->SetVertex(2, X_cm + width, Y_cm + height, Z_cm - length);
                            b->SetVertex(3, X_cm + width, Y_cm - height, Z_cm - length);
                            b->SetVertex(4, X_cm - width, Y_cm - height, Z_cm + length);
                            b->SetVertex(5, X_cm - width, Y_cm + height, Z_cm + length);
                            b->SetVertex(6, X_cm + width, Y_cm + height, Z_cm + length);
                            b->SetVertex(7, X_cm + width, Y_cm - height, Z_cm + length);
                        }
                        // Crv T orientation (Length along X)
                        else if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){
                            b->SetVertex(0, X_cm - length, Y_cm - height, Z_cm - width);
                            b->SetVertex(1, X_cm - length, Y_cm + height, Z_cm - width);
                            b->SetVertex(2, X_cm - length, Y_cm + height, Z_cm + width);
                            b->SetVertex(3, X_cm - length, Y_cm - height, Z_cm + width);
                            b->SetVertex(4, X_cm + length, Y_cm - height, Z_cm - width);
                            b->SetVertex(5, X_cm + length, Y_cm + height, Z_cm - width);
                            b->SetVertex(6, X_cm + length, Y_cm + height, Z_cm + width);
                            b->SetVertex(7, X_cm + length, Y_cm - height, Z_cm + width);
                        }
                        // Crv R, Crv L orientation (Length along Y)
                        else if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 1){
                            b->SetVertex(0, X_cm - height, Y_cm - length, Z_cm - width);
                            b->SetVertex(1, X_cm - height, Y_cm + length, Z_cm - width);
                            b->SetVertex(2, X_cm + height, Y_cm + length, Z_cm - width);
                            b->SetVertex(3, X_cm + height, Y_cm - length, Z_cm - width);
                            b->SetVertex(4, X_cm - height, Y_cm - length, Z_cm + width);
                            b->SetVertex(5, X_cm - height, Y_cm + length, Z_cm + width);
                            b->SetVertex(6, X_cm + height, Y_cm + length, Z_cm + width);
                            b->SetVertex(7, X_cm + height, Y_cm - length, Z_cm + width);
                        }
                        allcrvbars_collection->AddElement(b); 

                    } else if(extracted){
                        // T1 (Length along Z)
                        if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ 
                            b->SetVertex(0, X_cm - width, Y_cm - height, Z_cm - length);
                            b->SetVertex(1, X_cm - width, Y_cm + height, Z_cm - length);
                            b->SetVertex(2, X_cm + width, Y_cm + height, Z_cm - length);
                            b->SetVertex(3, X_cm + width, Y_cm - height, Z_cm - length);
                            b->SetVertex(4, X_cm - width, Y_cm - height, Z_cm + length);
                            b->SetVertex(5, X_cm - width, Y_cm + height, Z_cm + length);
                            b->SetVertex(6, X_cm + width, Y_cm + height, Z_cm + length);
                            b->SetVertex(7, X_cm + width, Y_cm - height, Z_cm + length);
                        } 
                        // EX, T2 (Length along X)
                        else { 
                            b->SetVertex(0, X_cm - length, Y_cm - height, Z_cm - width);
                            b->SetVertex(1, X_cm - length, Y_cm + height, Z_cm - width);
                            b->SetVertex(2, X_cm + length, Y_cm + height, Z_cm - width);
                            b->SetVertex(3, X_cm + length, Y_cm - height, Z_cm - width);
                            b->SetVertex(4, X_cm - length, Y_cm - height, Z_cm + width);
                            b->SetVertex(5, X_cm - length, Y_cm + height, Z_cm + width);
                            b->SetVertex(6, X_cm + length, Y_cm + height, Z_cm + width);
                            b->SetVertex(7, X_cm + length, Y_cm - height, Z_cm + width);
                        }
                        scene->AddElement(b); 
                    }
                }
            } // End of inner h loop (pulses)
            
            // Cluster Point Set Configuration
            ps1->SetMarkerColor(drawconfig.getInt("CrvHitColor")); // Marker color based on configuration
            ps1->SetMarkerStyle(DataInterface::mstyle);
            ps1->SetMarkerSize(DataInterface::msize);
            
            if(ps1->GetSize() !=0 ) scene->AddElement(ps1); 
            
        } // End of j loop (clusters)

        // Final Additions to Scene
        if(!extracted && addCrvBars) {
            scene->AddElement(allcrvbars_collection); 
        } 
    }
}

/*
 * Adds reconstructed TimeClusers data products to the REve visualization scene.
 * visualized as points
*/
// FIXME - do we ever use this? What is the purpose?
void DataInterface::AddTimeClusters(REX::REveManager *&eveMng, bool firstLoop_, 
                                 std::tuple<std::vector<std::string>, 
                                 std::vector<const TimeClusterCollection*>> timecluster_tuple, 
                                 REX::REveElement* &scene){

    
    std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(timecluster_tuple);
    std::vector<std::string> names = std::get<0>(timecluster_tuple);

    if(timecluster_list.empty()){
        return; 
    } else {
      std::cout << "[DataInterface::AddTimeClusters()]" << std::endl;
    }
    
    // Loop over TimeCluster Collections (i-loop)
    for(unsigned int i = 0; i < timecluster_list.size(); i++){
        const TimeClusterCollection* tccol = timecluster_list[i];
        
        if(tccol && tccol->size() != 0){
            
            // Create a compound for this entire collection for better organization in the REve tree
            std::string collection_title = "Time Cluster Collection: " + names[i];
            auto collection_compound = new REX::REveCompound(collection_title.c_str(), collection_title.c_str(), true);
            
            // Loop over individual TimeClusters (j-loop)
            for(size_t j=0; j<tccol->size();j++){
                mu2e::TimeCluster const &tclust= (*tccol)[j];
                
                // Detailed title for mouseover
                std::string tctitle = "Time Cluster tag: " + names[i] + '\n'
                + "t0: " + std::to_string(tclust.t0().t0()) + " +/- " + std::to_string(tclust.t0().t0Err()) + " ns " + '\n' 
                + "Hits: " + std::to_string(tclust.hits().size()) + '\n'
                + "Position (X, Y, Z): (" + std::to_string(tclust._pos.x()) + ", " + std::to_string(tclust._pos.y()) + ", " + std::to_string(tclust._pos.z()) + ") mm";
                
                // Create a point set for the single cluster. Naming convention improved.
                std::string ps_name = "Cluster_" + std::to_string(j) + "_" + names[i];
                auto ps1 = new REX::REvePointSet(ps_name.c_str(), tctitle.c_str(), 1); 
                
                // Get position and convert to cm
                CLHEP::Hep3Vector clusterPos(tclust._pos.x(), tclust._pos.y(), tclust._pos.z());
                ps1->SetNextPoint(pointmmTocm(clusterPos.x()), pointmmTocm(clusterPos.y()) , pointmmTocm(clusterPos.z()));
                
                // Marker style configuration
                ps1->SetMarkerColor(i + 6); // Use a color based on collection index (i)
                ps1->SetMarkerStyle(kOpenCircle);
                ps1->SetMarkerSize(DataInterface::msize * 1.5); // Slightly larger marker for clarity

                if(ps1->GetSize() !=0 ) {
                    // Add the individual cluster point set to the collection compound
                    collection_compound->AddElement(ps1);
                }
            } // End of j-loop (clusters)
            
            // Add the entire collection compound to the main scene
            scene->AddElement(collection_compound);
        }
    } // End of i-loop (collections)
}

/*
 * Adds reconstructed HelixSeed data products to the REve visualization scene.
 * Visualized as series of lines
*/
//FIXME - establish a use case, if none, remove!
void DataInterface::AddHelixSeedCollection(REX::REveManager *&eveMng, bool firstloop, 
                                         std::tuple<std::vector<std::string>, 
                                         std::vector<const HelixSeedCollection*>> helix_tuple, 
                                         REX::REveElement* &scene)
{
    std::cout << "[DataInterface::AddHelixSeedCollection]" << std::endl;

    //FIXME - remove hardcoded (OK, it is the size of the tracker)
    const int Z_START_MM = -1500;
    const int Z_END_MM   = 1500;
    const int Z_STEP_MM  = 100;
    const int NUM_POINTS = (Z_END_MM - Z_START_MM) / Z_STEP_MM + 1; 

    const auto& helix_list = std::get<1>(helix_tuple);
    const auto& names = std::get<0>(helix_tuple);
    
    if (helix_list.empty()) return;

    for(unsigned int j=0; j < helix_list.size(); j++){
        const HelixSeedCollection* seedcol = helix_list[j];
        
        if(seedcol != nullptr){
            std::string collection_title = "HelixSeed Collection: " + names[j];
            auto collection_compound = new REX::REveCompound(collection_title.c_str(), collection_title.c_str(), true);

            for(unsigned int k = 0; k < seedcol->size(); k++){
                const mu2e::HelixSeed& hseed = (*seedcol)[k];
                const mu2e::RobustHelix& rhelix = hseed.helix();
                
                // Calculate/Extract Useful Parameters for the Title
                double momentum = rhelix.momentum(); // Momentum magnitude
                double lambda_pitch = rhelix._lambda; // Pitch parameter: dZ/dPhi
                double r_center = rhelix._rcent; // Radius of center vector
                
                // Create Detailed Title String
                std::string helix_title = "HelixSeed " + std::to_string(k) + " tag: " + names[j] + '\n'
                    + " Momentum (p): " + std::to_string(momentum) + " MeV/c" + '\n'
                    + " Pitch (lambda): " + std::to_string(lambda_pitch) + '\n'
                    + " Center Radius (R_cent): " + std::to_string(r_center) + " mm";
                
                
                // Use the calculated number of points for the REveLine constructor
                auto line = new REX::REveLine(helix_title.c_str(), helix_title.c_str(), NUM_POINTS);
                
                // Get Helix Parameters for Drawing
                float helrad  = rhelix._radius;
                float xc      = rhelix._rcent * cos(rhelix._fcent);
                float yc      = rhelix._rcent * sin(rhelix._fcent);
                float lambda  = rhelix._lambda;
                float fz0     = rhelix._fz0;
                
                // Draw Synthetic Helix Points
                for(int i = Z_START_MM; i <= Z_END_MM; i += Z_STEP_MM){
                    float z = static_cast<float>(i);
                    float circphi = 0.0;
                    
                    if(lambda != 0.0f) {
                        circphi = fz0 + z / lambda;
                    } else {
                        circphi = fz0; 
                    }

                    float x = xc + helrad * cos(circphi);
                    float y = yc + helrad * sin(circphi);
                    
                    CLHEP::Hep3Vector HelPos(x, y, z);
                    
                    line->SetNextPoint(pointmmTocm(HelPos.x()), pointmmTocm(HelPos.y()), pointmmTocm(HelPos.z()));
                }
                
                // Update styles
                line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
                line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
                
                collection_compound->AddElement(line);
            } // End of k-loop (HelixSeeds)
            
            scene->AddElement(collection_compound);
        }
    }
}

/*
 * Adds reconstructed KalIntersections data products to the REve visualization scene.
 * Visualizedas points
*/
void DataInterface::AddKalIntersection(mu2e::KalSeed const& kalseed, REX::REveElement* &scene, REX::REveCompound *products, std::string track_tag)
{
    // Retrieve the vector of intersections from the KalSeed
    std::vector<mu2e::KalIntersection> const& inters = kalseed.intersections();

    // Iterate over all intersections
    for(mu2e::KalIntersection const& inter : inters){
        
        const KinKal::VEC3& posKI = inter.position3();
        
        // Create Detailed Title String
        std::string title = std::string("KalIntersection for track: ") + track_tag + '\n'
            + "Surface: " + inter.surfaceId().name() + '\n'
            + "Position (x, y, z): (" 
                + std::to_string(posKI.x()) + ", " 
                + std::to_string(posKI.y()) + ", " 
                + std::to_string(posKI.z()) + ") mm" + '\n'
            + "Time: " + std::to_string(inter.time()) + " ns" + '\n'
            + "Momentum (p, dp): " 
                + std::to_string(inter.mom()) + " +/- " 
                + std::to_string(inter.dMom()) + " MeV/c";
        
        // Create and Style REve Point Set
        auto interpoint = new REX::REvePointSet(title.c_str(), title.c_str(), 1);

        // Assign color based on momentum uncertainty (dMom > 0 for material intersections)
        Color_t marker_color = (fabs(inter.dMom()) > 0.0) ? kViolet : kYellow;
        
        interpoint->SetMarkerStyle(DataInterface::mstyle);
        interpoint->SetMarkerSize(DataInterface::msize);
        interpoint->SetMarkerColor(marker_color);

        // Set Point Position
        interpoint->SetNextPoint(pointmmTocm(posKI.x()), 
                                 pointmmTocm(posKI.y()), 
                                 pointmmTocm(posKI.z()));
        
        // Add to Products Compound
        if (interpoint->GetSize() != 0) {
            products->AddElement(interpoint);
        }
    }
}

/*
 * Adds reconstructed TrkStrawHits associated with KalSeed data products to the REve visualization scene.
 * Visualized as points
*/
template<class KTRAJc> 
void DataInterface::AddTrkStrawHit(mu2e::KalSeed const& kalseed, 
                                 REX::REveElement* &scene,  
                                 std::unique_ptr<KTRAJc> &lhptr, 
                                 REX::REveCompound *trackproducts)
{
    std::cout << "[DataInterface::AddTrkStrawHit]" << std::endl;
    
    // Setup and Data Extraction
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    std::vector<mu2e::TrkStrawHitSeed> const& hits = kalseed.hits();

    // Configuration Constants
    const double N_SIGMA = 2.0; // Length of the error line in sigma
    const std::string SIGMA_TITLE = "+/-" + std::to_string(N_SIGMA) + " sigma";
    
    // Loop Over Track-Straw Hits
    for(mu2e::TrkStrawHitSeed const& tshs : hits){
        
        // Get Geometry and State
        const mu2e::Straw& straw = tracker->straw(tshs.strawId());
        
        // Setup the hit state based on parameters in the seed
        mu2e::WireHitState whs(mu2e::WireHitState::State(tshs._ambig),
                               mu2e::StrawHitUpdaters::algorithm(tshs._algo),
                               tshs._kkshflag);
        
        bool active = whs.active();
        bool usedrift = whs.driftConstraint();

        // Calculate Position and Error Vectors
        if(active){ 
            
            // Start position: position on the wire at the reference POCA (RUpOS)
            auto tshspos = XYZVectorF(straw.wirePosition(tshs._rupos));
            double hit_error(0.0);
            
            // Find direction perpendicular to the wire and the track direction (drift direction)
            // track direction at POCA
            auto tdir = lhptr->direction(tshs._ptoca); 
            // wire direction
            auto wdir = XYZVectorF(straw.wireDirection(tshs._rupos)); 
            // drift direction is perpendicular to the plane formed by wire and track
            auto ddir = wdir.Cross(tdir).Unit() * whs.lrSign(); 

            if(usedrift){
                // Move position out by the drift distance along the signed drift direction
                tshspos += tshs._rdrift * ddir;
                hit_error = tshs._sderr; // Use signed drift error
            } else {
                hit_error = tshs._uderr; // Use unsigned drift error
            }
            
            // Calculate the endpoints for the error line segment
            auto end1 = tshspos + N_SIGMA * hit_error * ddir;
            auto end2 = tshspos - N_SIGMA * hit_error * ddir;
            
            // Create Title and Compound
            std::string title = std::string("TrkStrawHitSeed:") + '\n'
                + " Position (x,y,z): (" 
                    + std::to_string(tshspos.x()) + ", " 
                    + std::to_string(tshspos.y()) + ", " 
                    + std::to_string(tshspos.z()) + ") mm" + '\n'
                + " Time: " + std::to_string(tshs.time()) + " ns" + '\n'
                + " EnergyDep: " + std::to_string(tshs.energyDep()) + " MeV" + '\n'
                + " Error: " + SIGMA_TITLE;

            // Compound holds the point and the error line together
            auto point_with_error = new REX::REveCompound(title.c_str(), "TrkStrawHitSeed", 1);
            
            // Create Point Marker
            auto trkstrawpoint = new REX::REvePointSet(title.c_str(), title.c_str(), 1);
            
            // Marker configuration
            trkstrawpoint->SetMarkerStyle(DataInterface::mstyle);
            trkstrawpoint->SetMarkerSize(DataInterface::msize);
            
            // Color logic: Redraw color if drift constraint wasn't used
            Color_t base_color = drawconfig.getInt("TrkHitColor");
            if (!usedrift) {
                base_color = drawconfig.getInt("TrkNoHitColor");
            }
            trkstrawpoint->SetMarkerColor(base_color);
            
            // Set the position of the hit marker
            trkstrawpoint->SetNextPoint(pointmmTocm(tshspos.x()), 
                                       pointmmTocm(tshspos.y()), 
                                       pointmmTocm(tshspos.z()));
                                       
            //  Error Line
            auto line = new REX::REveLine(("TrkStrawHit Error " + SIGMA_TITLE).c_str(), SIGMA_TITLE.c_str(), 2);
            
            // Set line endpoints (converted to cm)
            line->SetNextPoint(pointmmTocm(end1.x()), pointmmTocm(end1.y()), pointmmTocm(end1.z()));
            line->SetNextPoint(pointmmTocm(end2.x()), pointmmTocm(end2.y()), pointmmTocm(end2.z()));
            
            // Line color is the same as the base hit color
            line->SetLineColor(base_color); 

            // Add Elements to Compound and Compound to Products
            point_with_error->AddElement(trkstrawpoint);
            point_with_error->AddElement(line);
            
            trackproducts->AddElement(point_with_error);
        }
    }
}

/*
 * Adds reconstructed TrkCaloHitSeed products, visualized as points
*/
void DataInterface::AddTrkCaloHit(mu2e::KalSeed const& kalseed, REX::REveElement* &scene)
{
    std::cout<<"[DataInterface::AddTrkCaloHit]"<<std::endl;
    // The TrkCaloHitSeed is a member of the KalSeed
    const mu2e::TrkCaloHitSeed& caloseed = kalseed.caloHit();
    
    // The caloseed contains an art::Ptr to the CaloCluster
    art::Ptr<mu2e::CaloCluster> cluster = caloseed.caloCluster(); 

    // Check if a valid CaloCluster is associated with the track
    if (cluster) {
        // Extract Cluster Information
        CLHEP::Hep3Vector clusterPos = cluster->cog3Vector();
        double energy = cluster->energyDep();
        double time = cluster->time();
        // Get Geometry Handles
        mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
        GeomHandle<DetectorSystem> det;
        // Create Title String
        std::string cluster_title = std::string("TrkCaloHit CaloCluster") + '\n'
                                  + " Position (x,y,z): (" 
                                  + std::to_string(clusterPos.x()) + ", " 
                                  + std::to_string(clusterPos.y()) + ", " 
                                  + std::to_string(clusterPos.z()) + ") mm" + '\n'
                                  + " Energy: " + std::to_string(energy) + " MeV" + '\n'
                                  + " Time: " + std::to_string(time) + " ns";

        CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(cluster->diskID(),clusterPos);
        CLHEP::Hep3Vector pointInMu2e = det->toMu2e(crystalPos);

        // Create REve Point Set
        auto ps1 = new REX::REvePointSet(cluster_title.c_str(), cluster_title.c_str(), 1); 

        // Set Positions, Color, and Size
        if(cluster->diskID() == 0)    
            ps1->SetNextPoint(pointmmTocm(clusterPos.x()), pointmmTocm(clusterPos.y()) , abs(pointmmTocm(pointInMu2e.z())));
        if(cluster->diskID() == 1)
            ps1->SetNextPoint(pointmmTocm(clusterPos.x()), pointmmTocm(clusterPos.y()) , abs(pointmmTocm(pointInMu2e.z())));


        // Styling
        ps1->SetMarkerColor(drawconfig.getInt("TrkHitColor")); 
        ps1->SetMarkerStyle(kFullDiamond);
        ps1->SetMarkerSize(DataInterface::msize * 3.0);

        // Add to Scene
        if (ps1->GetSize() != 0) {
            scene->AddElement(ps1);
        }
    }
}

/*
 * Adds reconstructed KTRAJ data products to the REve visualization scene.
 * Visualized as line set
*/
using LHPT = KalSeed::LHPT;
using CHPT = KalSeed::CHPT;
using KLPT = KalSeed::KLPT;
template<class KTRAJ> 
void DataInterface::AddKinKalTrajectory(std::unique_ptr<KTRAJ> &trajectory, 
                                        REX::REveElement* &scene, 
                                        unsigned int j, 
                                        std::string kaltitle, 
                                        double& t1, // event range will be updated
                                        double& t2)
{
    // Extract Time Range
    // The range object is returned by value/reference, so we access its limits once.
    t1 = trajectory->range().begin();
    t2 = trajectory->range().end();

    // Calculate Number of Points and Setup Line Object
    // Calculate required number of points for the loop: (t2 - t1) / 0.1 + 1 (for the starting point)
    double time_step = 0.1;
    size_t num_steps = static_cast<size_t>((t2 - t1) / time_step) + 1;

    // Create the REveLine with the pre-calculated title and required size
    auto line = new REX::REveLine(kaltitle.c_str(), kaltitle.c_str(), num_steps);

    // Iterate Through Time and Plot Points 
    
    // Set the first point explicitly (t = t1)
    const auto &p_start = trajectory->position3(t1);
    line->SetPoint(0, pointmmTocm(p_start.x()), pointmmTocm(p_start.y()), pointmmTocm(p_start.z()));
    
    // Loop from t1 + step up to t2
    for(double t = t1 + time_step; t <= t2; t += time_step)
    {
        // Get the position vector once
        const auto &p = trajectory->position3(t);
        
        // Add the point, converting units
        line->SetNextPoint(pointmmTocm(p.x()), 
                           pointmmTocm(p.y()), 
                           pointmmTocm(p.z()));
    }

    // Styling and Scene Addition
    line->SetLineColor(j + 6); // Use a color based on the collection index (j)
    line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
    scene->AddElement(line);
}

/*
 * Adds reconstructed KTRAJ
*/
void DataInterface::FillKinKalTrajectory(REX::REveManager *&eveMng, bool firstloop, REX::REveElement* &scene, 
                                         std::tuple<std::vector<std::string>, 
                                         std::vector<const KalSeedPtrCollection*>> track_tuple, 
                                         bool plotKalIntersection, bool addTrkHits, bool addTrkCaloHits, 
                                         double& t1, double& t2)
{
    std::cout << "[DataInterface::FillKinKalTrajectory()]" << std::endl;
    
    // Critical Logic: Scene Cleanup 
    if (!firstloop) {
        scene->DestroyElements();
    }

    // Setup and Data Extraction 
    const auto& ptable = GlobalConstantsHandle<ParticleDataList>();
    const auto& track_list = std::get<1>(track_tuple);
    const auto& names = std::get<0>(track_tuple);
    
    if (track_list.empty()) return;

    // Loop over KalSeed Collections 
    for(unsigned int j = 0; j < track_list.size(); j++){
        const mu2e::KalSeedPtrCollection* seedcol = track_list[j];
        
        if(seedcol == nullptr) continue;
        
        // Loop over individual KalSeeds 
        for(const auto& kseedptr : *seedcol){
            if (!kseedptr) continue;
            
            const mu2e::KalSeed& kseed = *kseedptr;
            std::string particle_name = ptable->particle(kseed.particle()).name();

            // Determine active hits for display info 
            unsigned nactive = 0;
            for (auto const& hit : kseed.hits()){ 
                if (hit.strawHitState() > mu2e::WireHitState::inactive) {
                    ++nactive; 
                }
            }
            
            // Find t0 and Setup Output Container
            double t0;
            kseed.t0Segment(t0);

            // Container for all points, hits, and intersections for THIS track
            std::string compound_name = "Track Products for KalSeed " + particle_name;
            REX::REveCompound *trackproducts = new REX::REveCompound(compound_name.c_str(), compound_name.c_str(), 1);
            
            std::stringstream ksstream;
            std::string tag = particle_name;

            // Loop Helix Fit 
            if(kseed.loopHelixFit()) {
                tag += " Loop Helix";
                auto trajectory = kseed.loopHelixFitTrajectory();
                const auto& lh = trajectory->nearestPiece(t0);
                auto momvec = lh.momentum3(t0);
                
                ksstream << particle_name << " LoopHelix "
                    << std::setw(6) << std::setprecision(3) << momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << '\n'
                    << "t0 " << lh.t0() << " ns, "
                    << "lam " << lh.lam() << " mm, "
                    << "rad " << lh.rad() << " mm" << '\n'
                    << "cx " << lh.cx() << " mm, "
                    << "cy " << lh.cy() << " mm, "
                    << "phi0 " << lh.phi0() << " rad" << '\n'
                    << "N active hits: " << nactive << ", Fit cons: " << kseed.fitConsistency() << '\n'
                    << "Instance: " << names[j];
                
                AddKinKalTrajectory<LHPT>(trajectory, scene, j, ksstream.str(), t1, t2);
                if(addTrkHits) {
                    AddTrkStrawHit<LHPT>(kseed, scene, trajectory, trackproducts);
                    if(addTrkCaloHits) AddTrkCaloHit(kseed, scene);
                }
            }
            
            // Central Helix Fit 
            else if(kseed.centralHelixFit()) {
                tag += " Central Helix";
                auto trajectory = kseed.centralHelixFitTrajectory();
                const auto& ch = trajectory->nearestPiece(t0);
                auto momvec = ch.momentum3(t0);
                
                ksstream << particle_name << " CentralHelix "
                    << std::setw(6) << std::setprecision(3) << momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << '\n'
                    << "t0 " << ch.t0() << " ns, "
                    << "tandip " << ch.tanDip() << '\n'
                    << "d0 " << ch.d0() << " mm, "
                    << "z0 " << ch.z0() << " mm, "
                    << "phi0 " << ch.phi0() << " rad" << '\n'
                    << "omega " << ch.omega() << " mm^-1" << '\n'
                    << "Track arrival time " << t1 << '\n'
                    << "Instance: " << names[j];
                
                AddKinKalTrajectory<CHPT>(trajectory, scene, j, ksstream.str(), t1, t2);
                if(addTrkHits) {
                    AddTrkStrawHit<CHPT>(kseed, scene, trajectory, trackproducts);
                    if(addTrkCaloHits) AddTrkCaloHit(kseed, scene);
                }
            }
            
            // Kinematic Line Fit 
            else if(kseed.kinematicLineFit()) {
                tag += " Kinematic Line";
                auto trajectory = kseed.kinematicLineFitTrajectory();
                const auto& kl = trajectory->nearestPiece(t0);
                auto momvec = kl.momentum3(t0);
                
                ksstream << particle_name << " KinematicLine "
                    << std::setw(6) << std::setprecision(3) << momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << '\n'
                    << "t0 " << kl.t0() << " ns, "
                    << "d0 " << kl.d0() << " mm, "
                    << "z0 " << kl.z0() << " mm" << '\n'
                    << "phi0 " << kl.phi0() << " rad, "
                    << "theta " << kl.theta() << " rad" << '\n'
                    << "Track arrival time " << t1 << '\n'
                    << "Instance: " << names[j];
                
                AddKinKalTrajectory<KLPT>(trajectory, scene, j, ksstream.str(), t1, t2);
                if(addTrkHits) {
                    AddTrkStrawHit<KLPT>(kseed, scene, trajectory, trackproducts);
                    if(addTrkCaloHits) AddTrkCaloHit(kseed, scene);
                }
            }
            
            // Final Additions for the Current KalSeed 
            if(plotKalIntersection) {
                AddKalIntersection(kseed, scene, trackproducts, tag);
            }
            
            // Add the compound of all hits/intersections if we added anything to it
            // The logic implies we add it if hits OR intersections were requested.
            if(addTrkHits || plotKalIntersection) {
                scene->AddElement(trackproducts);
            }

        } // End of inner loop (KalSeeds)
    } // End of outer loop (Collections)
}

/*
 * Adds reconstructed CosmicTrackSeed product, visualized as line
*/
/*Function to visualize CosmicTrackSeed fits (straight lines)-*/
void DataInterface::AddCosmicTrackFit(REX::REveManager *&eveMng, bool firstLoop_, 
                                      const mu2e::CosmicTrackSeedCollection *cosmiccol, 
                                      REX::REveElement* &scene)
{
    std::cout << "[DataInterface] AddCosmicTrackSeed" << std::endl;

    
    if (cosmiccol == nullptr) return;
    
    // Create a compound to hold all individual cosmic tracks for organization
    auto all_tracks_compound = new REX::REveCompound("Cosmic Tracks", "Cosmic Tracks", true);
    
    // Loop over individual CosmicTrackSeeds 
    for(size_t i = 0; i < cosmiccol->size(); i++){
        const mu2e::CosmicTrackSeed& sts = (*cosmiccol)[i];
        const mu2e::CosmicTrack& st = sts._track;
        
        // Ensure there are enough hits to define a track segment
        if (sts._straw_chits.size() < 2) continue;

        // Define Track Segment Endpoints 
        // The track is parameterized as: X(y) = A0 - A1*y, Z(y) = B0 - B1*y.
        // We plot the line segment between the Y-positions of the first and last hits.
        
        const mu2e::ComboHit& first_hit = sts._straw_chits.front();
        const mu2e::ComboHit& last_hit = sts._straw_chits.back();
        
        double ty1 = first_hit.pos().y();
        double ty2 = last_hit.pos().y();
        
        // Calculate the (X, Z) coordinates at Y1 and Y2 using the fit parameters
        double tx1 = st.MinuitParams.A0 - st.MinuitParams.A1 * ty1;
        double tx2 = st.MinuitParams.A0 - st.MinuitParams.A1 * ty2;
        double tz1 = st.MinuitParams.B0 - st.MinuitParams.B1 * ty1;
        double tz2 = st.MinuitParams.B0 - st.MinuitParams.B1 * ty2;
        
        // Create and Populate REveLine for the Current Track 
        std::string track_title = "Cosmic Track " + std::to_string(i) 
                                + " N_Hits: " + std::to_string(sts._straw_chits.size());
        auto line = new REX::REveLine(track_title.c_str(), track_title.c_str(), 2);
        
        // Set Point 1 (Start)
        line->SetNextPoint(pointmmTocm(tx1), pointmmTocm(ty1), pointmmTocm(tz1));
        
        // Set Point 2 (End)
        line->SetNextPoint(pointmmTocm(tx2), pointmmTocm(ty2), pointmmTocm(tz2));

        // 3. Styling and Addition 
        line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
        line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
        
        all_tracks_compound->AddElement(line);
    }
    
    // Add the compound of all tracks to the scene
    scene->AddElement(all_tracks_compound);
}

