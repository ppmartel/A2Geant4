/***** Time Projection Chamber
A Postuma 2021 *****/

#include "A2TPC.hh"

//then, include a bunch of useful Geant4 stuff
#include "G4SDManager.hh" //manage sensitive detector
#include "G4Tubs.hh" //cylinder
#include "G4Box.hh" //box
#include "G4Cons.hh" //cone
#include "G4Trd.hh" //trapezoid
#include "G4NistManager.hh"  //element manager
#include "G4VisAttributes.hh" //visualization
#include "G4PVPlacement.hh" //placement
#include "G4UnionSolid.hh" //union of several solide
#include "G4SubtractionSolid.hh" //subtract other solids
#include "G4RotationMatrix.hh" //rotate and place solids
#include "G4OpticalSurface.hh" //optical properties
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh" //units
#include "G4PhysicalConstants.hh" //constants
#include "A2SD.hh" //sensitive detector
#include "A2VisSD.hh" //visual sensitive detector
#include "A2MagneticField.hh" //magnetic fields
#include "A2ElectricField.hh" //my electric field class
//#include "A2EMField.hh"

#include "G4FieldManager.hh" //manage the fields
#include "G4TransportationManager.hh" //transport through fields
#include "CLHEP/Units/SystemOfUnits.h" //units
#include "G4IntegrationDriver.hh"

using namespace CLHEP;

//first, a constructor
A2TPC::A2TPC(){
	//initiate all pointers to NULL
	
	//first, volumes: the three that are in every tgt file
	fMotherLogic = NULL; //mother logical volume
	fMyLogic = NULL; //the constructed logical volume
	fMyPhysi = NULL; //the constructed physical volume
	
	//for this detector class
	//volumes in MakeVessel()
	fVesselLogic=NULL;
	fMainCellLogic = NULL;
	//fMainCellHeLogic = NULL;
	fMainCellEndLogic = NULL;
	fExtCellLogic = NULL;
	//fExtCellHeLogic = NULL;
	fConeLogic = NULL;
	//fConeHeLogic = NULL;
	fBeWindowLogic = NULL;
	fVesselHeLogic = NULL;

	//volumes in MakeAnodeCathode()
	fAnodeLogic = NULL;
        fAnodeCentreLogic= NULL;
        fAnodeRingLogic= NULL;
	fCathodeLogic= NULL;
	for(G4int i=0;i<6;i++){
		fAnodeSecLogic[i]=NULL;
	}
	
	//volumes in MakeGrid()
	fGridLogic= NULL;
	for(G4int j=0;j<200;j++){
		fWireLogic[j]=NULL;
	}
	
	//some sensitive detector stuff
        fRegionAnode = new G4Region("Anode");
        fAnodeSD = NULL;
        fAnodeVisSD = NULL;

	//some E field stuff
	fElectricField = NULL;
	
	//initiate the nist manager
	fNistManager=G4NistManager::Instance();

	//initiate other relevant pointers here
	fIsOverlapVol=true;

	//read dimensions from a parameter file
	ReadParameters("data/TPC.dat");
}

//destructor
A2TPC::~A2TPC(){
	if(fRegionAnode) delete fRegionAnode; //remove the sensitive detector
	if(fAnodeSD) delete fAnodeSD; //remove the sensitive detector
	if(fAnodeVisSD) delete fAnodeVisSD; //remove the sensitive detector
	if(fElectricField) delete fElectricField; //remove electric field
}

//main constructor function called inside main construction
//maybe where the pointer is initiated to NULL??
G4VPhysicalVolume* A2TPC::Construct(G4LogicalVolume* MotherLogical, G4double Z0){
	
	//get pointer to mother logical volume
	fMotherLogic=MotherLogical;
	
	//output a message
	G4cout<< "A2TPC::Construct() Building the TPC." <<G4endl;
	
	//call functions
	DefineMaterials(); //define relevant materials
	MakeVessel(); //build fVesselLogic
	MakeAnodeCathode(); //build fAnodeLogic and fCathodeLogic
	MakeGrid();
	MakeSensitiveDetector(); //make the anode into a detector
	MakeField(); //create field inside the detector
	PlaceParts(); //place the seperate detector parts into fMyLogic
	
	//place fMyLogic into fMotherLogic
	fMyPhysi = new G4PVPlacement(0, G4ThreeVector(0,0,Z0) ,fMyLogic,"TPC",fMotherLogic,false,1,fIsOverlapVol);

	//and return this physical volume describing the target!
	return fMyPhysi;
}

/***** This function builds the target vessel and fills the target with helium *****/
void A2TPC::MakeVessel(){
	/***** shapes ******/
	//overall vessel shape
	G4Tubs *fVessel = new G4Tubs("Vessel",
				0,
				fRadius+fThickness,
				fLength/2+fExtension+fConeLength,
				0.*deg,
				360.*deg);
	//main cell steel shell
	G4Tubs *fMainCell = new G4Tubs("MainCell", //name
				fRadius,
				fRadius+fThickness,
				fLength/2,
                                0.*deg,       //start angle
                                360.*deg);    //spanning angle
	//circular end of main cell, with hole for extension
	G4Tubs *fMainCellEnd = new G4Tubs("MainCellEnd",
				fExtRadius+fThickness,
				fRadius+fThickness,
				fEndThickness/2,
				0.*deg,
				360.*deg);
	//conical main cell end on the other side
        G4Cons *fMainCellCone = new G4Cons("MainCellCone",
                                fRadius, //inner radius at one end
                                fRadius+fThickness, //outer
                                fExtRadius, //inner radius at other end
                                fExtRadius+fThickness, //outer
                                fConeLength/2, //half length
                                0.*deg, //start angle
                                360.*deg); //spanning angle
	//helium inside main cell
	G4Tubs *fMainCellHe = new G4Tubs("MainCellHelium",
				0,
				fRadius,
				fLength/2,
				0.*deg,
				360.*deg);
        //helium inside conical part of cell
	G4Cons *fConeHe = new G4Cons("ConeHe",
                                0,
                                fRadius,
                                0,
                                fExtRadius,
                                fConeLength/2,
                                0.*deg,
                                360.*deg);
	//extension tube
	//made once and placed twice
	G4Tubs *fExtCell = new G4Tubs("ExtCell",
				fExtRadius,
				fExtRadius+fThickness,
				fExtension/2,
				0.*deg,
				360.*deg);
	//helium inside main extension tubes
	G4Tubs *fExtCellHe1 = new G4Tubs("ExtCellHe",
				0,
				fExtRadius,
				fExtension/2,
				0.*deg,
				360.*deg);
	G4Tubs *fExtCellHe2 = new G4Tubs("ExtCellHe",
				0,
				fExtRadius,
				fExtension/2+fEndThickness/2,
				0.*deg,
				360.*deg);
	//beryllium window at the end of each extension tube
        G4Tubs *fBeWindow = new G4Tubs("BeWindow",
                                0,
                                41.*mm,
                                fBeThickness/2,
                                0.*deg,
                                360.*deg);
	
	/***** create union solid of all helium *****/	
	G4UnionSolid *fHeUnion1 = new G4UnionSolid("HeUnion1",
			fMainCellHe,
			fConeHe,
			0,
			G4ThreeVector(0,0,(fLength+fConeLength)/2));
	G4UnionSolid *fHeUnion2 = new G4UnionSolid("HeUnion2",
			fHeUnion1,
			fExtCellHe1,
			0,
			G4ThreeVector(0,0,(fLength+fExtension)/2+fConeLength));
	G4UnionSolid *fVesselHe = new G4UnionSolid("VesselHe",
			fHeUnion2,
			fExtCellHe2,
			0,
			G4ThreeVector(0,0,-(fLength+fExtension+fEndThickness)/2));
	
	/**** logical volumes ******/
	//vessel
	fVesselLogic = new G4LogicalVolume
		(fVessel,
		 fNistManager->FindOrBuildMaterial("G4_AIR"),
		 "VesselLogic");
	//main cell
	fMainCellLogic = new G4LogicalVolume
            	(fMainCell,   //solid
            	 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"), //material
             "MainCellLogic"); //name for logical volume
	//main cell endpiece
	fMainCellEndLogic = new G4LogicalVolume
		(fMainCellEnd,
		 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
		 "MainCellEndLogic");
	//conical end to main cell
        fConeLogic = new G4LogicalVolume
                (fMainCellCone,
                 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
                 "ConeLogic");
	//helium
	fVesselHeLogic = new G4LogicalVolume
		(fVesselHe,
		 fNistManager->FindOrBuildMaterial("ATGasPure"),
		 "VesselHeLogic");
	//extension cells
	fExtCellLogic = new G4LogicalVolume
		(fExtCell,
		 fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
		 "ExtCellLogic");
	//beryllium windows
	fBeWindowLogic = new G4LogicalVolume
                (fBeWindow,
                 fNistManager->FindOrBuildMaterial("ATBerylliumW"),
                 "BeWindowLogic");
	

	/**** set visualization attributes *****/
	G4VisAttributes* lblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
    	G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
    	G4VisAttributes* cyan   = new G4VisAttributes( G4Colour(0.0,1.0,1.0,0.5)  );
	fVesselLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fMainCellLogic->SetVisAttributes(grey);
	//fMainCellLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fVesselHeLogic->SetVisAttributes(cyan);
	//fVesselHeLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fMainCellEndLogic->SetVisAttributes(grey);
    	fConeLogic->SetVisAttributes(grey);
	//fMainCellEndLogic->SetVisAttributes(G4VisAttributes::Invisible);
    	//fConeLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fExtCellLogic->SetVisAttributes(grey);
	//fExtCellLogic->SetVisAttributes(G4VisAttributes::Invisible);
	fBeWindowLogic->SetVisAttributes(lblue);
	//fBeWindowLogic->SetVisAttributes(G4VisAttributes::Invisible);

	/**** place things inside of VesselLogic ****/
	//main cell
    	new G4PVPlacement (0,                    //rotation
                     G4ThreeVector(0,0,0), //position
                     fMainCellLogic,         //logical volume
                     "MainCellPlacement",       //name
                     fVesselLogic,             //mother logic
                     false,                //pMany = false always
                     2,fIsOverlapVol);    //unique copy number
	//main cell endpiece
    	new G4PVPlacement(0,
		    G4ThreeVector(0,0,-(fLength+fEndThickness)/2),
		    fMainCellEndLogic,
		    "MainCellEndPlacement",
		    fVesselLogic,
		    3,fIsOverlapVol);
    	//main cell cone
	new G4PVPlacement(0,
		    G4ThreeVector(0,0,(fLength+fConeLength)/2),
		    fConeLogic,
		    "ConePlacement",
		    fVesselLogic,
		    4,fIsOverlapVol);
	//extension cells
	new G4PVPlacement (0,
		   	G4ThreeVector(0,0,(fLength+fExtension)/2+fConeLength),
			fExtCellLogic,
			"ExtCellPlacement1",
			fVesselLogic,
			6,fIsOverlapVol);
	new G4PVPlacement(0,
			G4ThreeVector(0,0,-(fLength+fExtension)/2-fEndThickness),
			fExtCellLogic,
			"ExtCellPlacement2",
			fVesselLogic,
			7,fIsOverlapVol);
	//beryllium windows - for now overlapping with helium
        new G4PVPlacement(0,
                        G4ThreeVector(0,0,fLength/2 + fExtension -fBeThickness/2+fConeLength),
                        fBeWindowLogic,
                        "BeWindowPlacement1",
                        fVesselLogic,
                        10,fIsOverlapVol);
        new G4PVPlacement(0,
                        G4ThreeVector(0,0,-(fLength/2+fExtension-fBeThickness/2+fEndThickness)),
                        fBeWindowLogic,
                        "BeWindowPlacement2",
                        fVesselLogic,
                        11,fIsOverlapVol);
}

/***** This function creates a segmented anode and a cathode *****/
void A2TPC::MakeAnodeCathode(){
        /***** solids for anode geometry *****/
        //main volume to hold sub-pieces
	G4Tubs *fAnode = new G4Tubs("Anode",
                                        0,
                                        fRadius,
                                        (fGThickness)/2,
                                        0.*deg,
                                        360.*deg);
	//circular central piece (G-10)
        G4Tubs *fAnodeCentre = new G4Tubs("AnodeCentre",
                                        0,
                                        10.*mm,
                                        fGThickness/2,
                                        0.*deg,
                                        360.*deg);
	//ring around central piece (G-10)
        G4Tubs *fAnodeRing = new G4Tubs("AnodeRing", 
					10.*mm,
                                        20.*mm,
                                        fGThickness/2,
                                        0.*deg,
                                        360.*deg);
	
	/***** logical volumes *****/	
	//volume holding entire anode
	fAnodeLogic = new G4LogicalVolume
                        (fAnode,
                         fNistManager->FindOrBuildMaterial("ATGasPure"), //this one contains the others - should be helium like the rest of the cell
                         "AnodeLogic");
	//circular central piece (G-10)
        fAnodeCentreLogic = new G4LogicalVolume
                        (fAnodeCentre,
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeCentreLogic");
	//ring around central piece (G-10)
        fAnodeRingLogic = new G4LogicalVolume
                        (fAnodeRing,
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeRingLogic");
	
	/***** visualization attributes *****/
        G4VisAttributes* lblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
        G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
        fAnodeLogic->SetVisAttributes(G4VisAttributes::Invisible);
        fAnodeCentreLogic->SetVisAttributes(lblue);
        fAnodeRingLogic->SetVisAttributes(lblue);

        /***** place anode parts to fAnodeLogic *****/
        //circular central piece (G-10)
	new G4PVPlacement(0, //rotation
                        G4ThreeVector(0,0,0), //placement
                        fAnodeCentreLogic, //logical volume
                        "AnodeCentrePlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        66,fIsOverlapVol); //unique copy number
        //ring around central piece (G-10)
	new G4PVPlacement(0, //rotation
                        G4ThreeVector(0,0,0), //placement
                        fAnodeRingLogic, //logical volume
                        "AnodeRingPlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        65,fIsOverlapVol); //unique copy number
	
	/*** segments of anode: reads specifications from parameter file ****/
	/***** define solid and logical volume for each radial section *****/
	for(G4int k=0; k<fRadialSecs; k++){
		//section of first row (G-10)
		fAnodeSec[k] = new G4Tubs("AnodeSec1",
                                        20.*mm+(k)*(fAnodeRadius-20.*mm)/fRadialSecs,
                                        20.*mm+(k+1)*(fAnodeRadius-20.*mm)/fRadialSecs,
                                        fGThickness/2,
                                        0.*deg,
					360.*deg/fAngularSecs);
        	fAnodeSecLogic[k] = new G4LogicalVolume
                        (fAnodeSec[k],
                        fNistManager->FindOrBuildMaterial("G-10"),
                        "AnodeSecLogic");
		fAnodeSecLogic[k]->SetVisAttributes(lblue);

		/***** place ring segments to anode *****/
        	for(G4int h=0; h<fAngularSecs;h++){ //for however many angular sections were defined in paramter file
          		//define the angle
                	G4double theta = h*360.*deg/fAngularSecs;
                	//generate a rotation matrix
                	G4RotationMatrix* rot = new G4RotationMatrix(theta,0,0);
                	//place one section of each ring at this angle
			//section of first row (G-10)
                	new G4PVPlacement(rot, //rotation
                        G4ThreeVector(0,0,0),
                        fAnodeSecLogic[k], //logical volume
                        "AnodeSecPlacement", //name
                        fAnodeLogic, //mother volume
                        false, //pmany: always false
                        1+h+k*fAngularSecs,fIsOverlapVol); //unique copy number
		}
	}
        
	/***** build solids for cathode geometry *****/
        //aluminum cathode
	G4Tubs *fCathode = new G4Tubs("Cathode",
                                        0,
                                        fCathodeRadius,
                                        fAlThickness/2,
                                        0.*deg,
                                        360.*deg);
        
	/***** cathode logical volumes *****/
        //aluminum cathode
	fCathodeLogic = new G4LogicalVolume(fCathode,
                        fNistManager->FindOrBuildMaterial("Aluminum"),
                        "CathodeLogic");
        
	/***** visulization attributes *****/
	fCathodeLogic->SetVisAttributes(lblue);
}

/**** This function builds the wire grid *****/
void A2TPC::MakeGrid(){
	//main volume to contain wires
	//solid: cylinder
	G4Tubs* fGrid = new G4Tubs("Grid", //name
				0, //inner radius
				fRadius, //outer radius
				2.*mm, //half length in Z
				0, //start angle
				360.*deg); //spanning angle
	//logical volume
	fGridLogic = new G4LogicalVolume(
			fGrid, //solid
			fNistManager->FindOrBuildMaterial("ATGasPure"), //material: helium, same as fVesselHeLogic
			"GridLogic"); //name
	G4int nWires = fRadius/fWireSpacing-1; //how many different length wires to build
	G4double fHWL[200]; //half length of each wire
	G4VisAttributes* grey   = new G4VisAttributes( G4Colour(0.5,0.5,0.5)  );
        fGridLogic->SetVisAttributes(G4VisAttributes::Invisible);
	
	//make and place each individual wire
	for(G4int l=0; l<nWires; l++){
		//find length needed to cross cylinder at a certain dinstace from the origin
		fHWL[l] = sqrt(fRadius*fRadius - l*l*fWireSpacing*fWireSpacing);
		//define wire with correct length
		fWire[l] = new G4Tubs("Wire",
			0, //inner radius
			fWireThickness, //outer radius
			fHWL[l], //half length
			0.*deg, //start angle
			360.*deg); //spanning angle
		//make the logical volume
		fWireLogic[l] = new G4LogicalVolume(fWire[l],
			fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL"),
			"WireLogic");
        	fWireLogic[l]->SetVisAttributes(grey);
		//determine distance from origin
		G4double XYPos = l*fWireSpacing;
		//create rotation matrices to align wire with X or Y axis
		G4RotationMatrix *xrot = new G4RotationMatrix(0,90.*deg,0);
		G4RotationMatrix *yrot = new G4RotationMatrix(90.*deg,90.*deg,0);
		//place 4 copies of each wire:
		//pointing along Y, displaced from origin in +X
		new G4PVPlacement(
				xrot, //rotation
				G4ThreeVector(XYPos,0,fWireThickness), //spatial coordinates 
				fWireLogic[l], //logical volume
				"WirePlacementX", //name
				fGridLogic, //mother logical
				false, //pMany: always false
				143, //unique copy no
			       	false); //check overlaps: not here
		//pointing along Y, displaced from origin in +X
		new G4PVPlacement(
				xrot,
				G4ThreeVector(-XYPos,0,fWireThickness),
				fWireLogic[l],
				"WirePlacementX2",
				fGridLogic,
				false,
				143, false);
		//pointing along X, displaced from origin in +Y
		new G4PVPlacement(
				yrot,
				G4ThreeVector(0,XYPos,-fWireThickness),
				fWireLogic[l],
				"WirePlacementY",
				fGridLogic,
				false,
				144,false);
		//pointing along X, displaced from origin in -Y
		new G4PVPlacement(
				yrot,
				G4ThreeVector(0,-XYPos,-fWireThickness),
				fWireLogic[l],
				"WirePlacementY2",
				fGridLogic,
				false,
				144,false);
	}

}

/***** this function makes the anode into a sensitive detector *****/
void A2TPC::MakeSensitiveDetector(){
        if(!fAnodeSD){ //check if SD is already defined
        G4SDManager* SDman = G4SDManager::GetSDMpointer(); //get pointer to SD manager
        fAnodeSD = new A2SD("AnodeSD",2+fRadialSecs*fAngularSecs); //create a new sensitive detector
        SDman->AddNewDetector(fAnodeSD); //add this detector to the SD manager
        //set each piece of anode as part of the sensitive detector
	fAnodeCentreLogic->SetSensitiveDetector(fAnodeSD);
        fAnodeRingLogic->SetSensitiveDetector(fAnodeSD);
	//set each piece of anode as part of anode region
        fRegionAnode->AddRootLogicalVolume(fAnodeCentreLogic);
	fRegionAnode->AddRootLogicalVolume(fAnodeRingLogic);
	for(G4int n=0; n<fRadialSecs; n++){
		fAnodeSecLogic[n]->SetSensitiveDetector(fAnodeSD);
		fRegionAnode->AddRootLogicalVolume(fAnodeSecLogic[n]);
	}
    }
}

/***** this function creates a uniform electric field through the target *****/
void A2TPC::MakeField(){
	//this is getting unruly 
	fElectricField = new A2ElectricField(); //create a field at zero first
	fElectricField->Construct(2);
	//fElectricField->SetDetectorField(fVesselLogic);
	/***
	G4cout<<"Making TPC electric field..."<<G4endl;
	//create the field: for now uniform, strength 100 kV/cm
	//NTS: what direction: probably in z: plus or minus??? should it point towards or away from cathode???
	fElectricField = new G4UniformElectricField(G4ThreeVector(0.0,0.0,2.0*kilovolt/cm));
	//fElectricField = new G4UniformElectricField(G4ThreeVector(0.0,0.0,100.0*kilovolt/cm));
	fElectricField = new A2ElectricField(); //create a field at zero first
	//G4ThreeVector fieldVector(0.0,0.0,2*kilovolt/cm); //define a value to give it
	//fElectricField->SetPureElectricFieldValue(fieldVector); //assign field to be electric at this strength
	//G4FieldManager* fieldManager = fElectricField->GetGlobalFieldManager();
	//fieldManager->SetDetectorField(fElectricField);
	//get field manager
	//G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
	G4FieldManager* fieldManager = new G4FieldManager(fElectricField);
        fVesselLogic->SetFieldManager(fieldManager,true);
	//define equation of motion of field
	fEquation = new G4EqMagElectricField(fElectricField);
	//now a stepper: according to default Geant4 method of calculating stepper in E or B field
	fEquation->SetFieldObj(fElectricField);
	auto fStepper = new G4DormandPrince745(fEquation,8); //eq of motion, nvariables
	//and an integration driver, to perform the integration of the steps
	fIntegrationDriver = new G4IntegrationDriver<G4DormandPrince745>(0.010*mm, //min step
							fStepper, //stepper
							8); //nvariables
	//attach to chord finder: object that calculates the trajectory
	fChordFinder = new G4ChordFinder(fIntegrationDriver);
	fieldManager->SetChordFinder(fChordFinder);
	fieldManager->SetDetectorField(fElectricField);
	//manually set step size: should be a fraction of the physics step size
	//fieldManager->GetChordFinder()->SetDeltaChord( 0.5*mm ); // Units: length
	                          // Relative accuracy values:
	//G4double minEps= 1.0e-7;  //   Minimum & value for largest steps
	//G4double maxEps= 1.0e-6;  //   Maximum & value for smallest steps

	//fieldManager->SetMinimumEpsilonStep( minEps );
	//fieldManager->SetMaximumEpsilonStep( maxEps );
	//fieldManager->SetDeltaOneStep( 0.5e-5* mm );  // 0.5 micrometer
**/

}


/***** this function places all parts of the target to fMyLogic *****/
void A2TPC::PlaceParts(){
	/***** define fMyLogic *****/
	//create a simple cylinder with adequate room for main and extension cells
	G4Tubs *fMyLogicCyl = new G4Tubs
            ("MyLogicCyl",
             0,             //inner radius
             fRadius+fThickness,   //outer radius
             fLength/2+fExtension+25.*mm,  //length/2
             0.*deg,        //start angle
             360.*deg);     //final angle
	//turn it into a logical volume
	fMyLogic = new G4LogicalVolume
            (fMyLogicCyl,  //solid
             fNistManager->FindOrBuildMaterial("G4_AIR"), //material
             "MyLogic");    //name
	//and set it invisible
    	fMyLogic->SetVisAttributes(G4VisAttributes::Invisible);
    
	/**** place anode and cathode to fVesselHeLogic *****/ 
        //first anode, close to beam entrance
	new G4PVPlacement(0,
                        G4ThreeVector(0,0,-fLength/2+fAnodeDistance),
                        fAnodeLogic,
                        "ANODE",
                        fVesselHeLogic,
                        30,fIsOverlapVol);
	//then cathode, close to beam exit
        new G4PVPlacement(0,
                        G4ThreeVector(0,0,fLength/2-fCathodeDistance),
                        fCathodeLogic,
                        "CATHODE",
                        fVesselHeLogic,
                        40,fIsOverlapVol); 
	/**** place grid to fVesselHeLogic ****/ 
	new G4PVPlacement(0,
			G4ThreeVector(0,0,-fLength/2+fAnodeDistance+3.*mm),
			fGridLogic,
			"GRID",
			fVesselHeLogic,
			50,fIsOverlapVol);
	/**** place helium containing anode, cathode, and grid to fVesselLogic ****/
	new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			fVesselHeLogic,
			"HELIUM",
			fVesselLogic,
			9,fIsOverlapVol);
	/**** place vessel inside fMyLogic *****/
	new G4PVPlacement(0,
			G4ThreeVector(0,0,0),
			fVesselLogic,
			"VESSEL",
			fMyLogic,
			10,fIsOverlapVol);

}

/***** this function reads target dimensions from a data file *****/
void A2TPC::ReadParameters(const char* file){
	//set keys contained in file
        char* keylist[] = { (char*) "TPC-Dim:", (char*) "Anode-Dim:", (char*) "Cathode-Dim:", (char*) "Grid-Dim:", (char*) "Run-Mode:", NULL};
        enum { ETPC_dim, EAnode_dim, ECathode_dim, EGrid_dim, ERun_Mode, ENULL };
        //define variables needed for reading the file
        G4int ikey, iread;
        G4int ierr = 0;
        char line[256];
        char delim[64];
        //open the file
        FILE* pdata;
        if( (pdata = fopen(file,"r")) == NULL ){
                printf("Error opening detector parameter file: %s\n",file);
                return;
        }
        //start reading the file
        while( fgets(line,256,pdata) ){
        if( line[0] == '#' ) continue; //skip commented lines
        printf("%s\n",line);
        sscanf(line,"%s",delim);
        for(ikey=0; ikey<ENULL; ikey++) //look for the defined keys
            if(!strcmp(keylist[ikey],delim)) break;
        switch( ikey ){
        default:
		//stop running if the file contains something strange
            printf("Unrecognised delimiter: %s\n",delim);
	    ierr++;
            break;
        case ETPC_dim: //dimensions of main target vessel
            iread = sscanf(line,"%*s%lf%lf%lf%lf%lf%lf%lf%lf", //format of line
                           &fLength,&fRadius,&fThickness,&fConeLength, //pointers to associated variables
                           &fExtension,&fExtRadius,&fEndThickness,&fBeThickness);
            if( iread != 8) ierr++;
            break;
        case EAnode_dim: //dimensions of anode
            iread = sscanf(line,"%*s%lf%lf%lf%i%i",
                            &fGThickness,&fAnodeRadius,&fAnodeDistance,&fRadialSecs,&fAngularSecs);
            if (iread !=5) ierr++;
            break;
        case EGrid_dim: //dimensions of anode
            iread = sscanf(line,"%*s%lf%lf",
                            &fWireThickness,&fWireSpacing);
            if (iread !=2) ierr++;
            break;
        case ECathode_dim: //dimensions of cathode
            iread = sscanf(line,"%*s%lf%lf%lf",
                            &fAlThickness,&fCathodeRadius,&fCathodeDistance);
            if (iread !=3) ierr++;
            break;
        case ERun_Mode: //run mode
            iread = sscanf(line,"%*s%d",
                           &fIsOverlapVol);
            if( iread != 1 ) ierr++;
            break;
	    }
        if( ierr ){ //if something strange is in the file, error and exit
    		printf("Fatal Error: invalid read of parameter line %s\n %s\n",
                   keylist[ikey],line);
            exit(-1);
                }
        }
}


/***** this function defines the materials used to build the target ******/
void A2TPC::DefineMaterials()
{
	G4double density, fractionmass;
	G4int ncomponents;
	//G4double pressure, temperature, a, z;

	/***** beryllium used in target windows *****/
	G4Material* BerylliumW = new G4Material("ATBerylliumW", density = 1.8480*g/cm3, ncomponents = 7);
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(14), 0.06*perCent);      //Si
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(4), 98.73*perCent);      //Be
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(6), 0.15*perCent);       //C
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(8), 0.75*perCent);       //O
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(12), 0.08*perCent);      //Mg
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(13), 0.1*perCent);       //Al
	BerylliumW->AddElement(fNistManager->FindOrBuildElement(26), 0.13*perCent);      //Fe

	/***** anode and cathode materials ****/
        //copied from online source
        //Anode G-10
        G4Material* G10 = new G4Material("G-10", density= 1.700*g/cm3, ncomponents=4);
        G10->AddElement(fNistManager->FindOrBuildElement(14), 1); //silicon
        G10->AddElement(fNistManager->FindOrBuildElement(8), 2); //oxygen
        G10->AddElement(fNistManager->FindOrBuildElement(6), 3); //carbon
        G10->AddElement(fNistManager->FindOrBuildElement(1), 3); //hydrogen
        //Anode copper
        G4Material* Cu = new G4Material("Copper",density= 8.960*g/cm3,ncomponents=1);
        Cu->AddElement(fNistManager->FindOrBuildElement(29),1);
        //Cathode aluminum
        G4Material* Al = new G4Material("Aluminum", density= 2.700*g/cm3,ncomponents=1);
        Al->AddElement(fNistManager->FindOrBuildElement(13),1);
    
	
	/***** other stuff from a2ActiveHe3.cc *****/
	G4Material* MylarW = new G4Material("ATMylarW", density = 1.4000*g/cm3, ncomponents = 3);
	MylarW->AddElement(fNistManager->FindOrBuildElement(1), 4.2*perCent);           //H
	MylarW->AddElement(fNistManager->FindOrBuildElement(6), 62.5*perCent);          //C
	MylarW->AddElement(fNistManager->FindOrBuildElement(8), 33.3*perCent);          //O

	G4Material* Teflon = new G4Material("ATTeflon", density = 2.2000*g/cm3, ncomponents = 2);
	Teflon->AddElement(fNistManager->FindOrBuildElement(6), 24*perCent);            //C
	Teflon->AddElement(fNistManager->FindOrBuildElement(9), 76*perCent);            //F


	//--------------------------------------------------
	// WLSfiber PMMA
	//--------------------------------------------------

	G4Material* PMMA = new G4Material("PMMA", density = 1.190*g/cm3, ncomponents = 3);
	PMMA->AddElement(fNistManager->FindOrBuildElement(6), 5); //C
	PMMA->AddElement(fNistManager->FindOrBuildElement(1), 8); //H
	PMMA->AddElement(fNistManager->FindOrBuildElement(8), 2); //O

	//Gas mixture and He3 management----------------------------------------------------

	//defining He3 according to
	// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/hadronprocess/731/1.html
	G4Element* ATHe3 = new G4Element("ATHe3", "ATHe3", ncomponents=1);
	ATHe3->AddIsotope((G4Isotope*)fNistManager->FindOrBuildElement(2)->GetIsotope(0), //isot. 0 = 3he
                      100.*perCent);

	//See HeGasDensity.nb. Might be an idea to include N2 effect to density
	//G4double he3density = 0.00247621*g/cm3; //20bar, calculated from ideal gas law
	//G4double he3density = 0.0033*g/cm3; //20bar, calculated from ideal gas law
	//G4double he3density = 0.004125*g/cm3; //25bar, calculated from ideal gas law
	G4double he3density = 0.00495*g/cm3; //30bar, calculated from ideal gas law

	G4Material* GasMix = new G4Material("ATGasMix", he3density, ncomponents = 2);
	GasMix->AddElement(ATHe3, 99.95*perCent);                                       //He3
	// GasMix->AddElement(fNistManager->FindOrBuildElement(2), 99.95*perCent);         //He4
	GasMix->AddElement(fNistManager->FindOrBuildElement(7), 0.05*perCent);           //N
	//Gas mixture and He3 end-----------------------------------------------------
	//----! IMPORTANT! Epoxy CURRENTLY TAKEN FROM A2 SIMULATION,------------------
	//NO IDEA WHETHER IT IS CORRECT OR NOT

	G4Material* GasPure = new G4Material("ATGasPure", he3density, ncomponents = 1);
	GasPure->AddElement(ATHe3, 100.*perCent);                                       //He3

	//Epoxy resin (C21H25Cl05) ***not certain if correct chemical formula or density
	//for colder temperature***:
	G4Material* EpoxyResin=new G4Material("EpoxyResin", density=1.15*g/cm3, ncomponents=4);
	EpoxyResin->AddElement(fNistManager->FindOrBuildElement(6), 21); //carbon
	EpoxyResin->AddElement(fNistManager->FindOrBuildElement(1), 25); //hydrogen
	EpoxyResin->AddElement(fNistManager->FindOrBuildElement(17), 1); //chlorine
	EpoxyResin->AddElement(fNistManager->FindOrBuildElement(8), 5);  //oxygen

	//Amine Hardener (C8H18N2) ***not certain if correct chemical formula or density
	//for colder temperature***:
	G4Material* Epoxy13BAC=new G4Material("Epoxy13BAC", density=0.94*g/cm3, ncomponents=3);
	Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(6), 8);
	Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(1), 18);
	Epoxy13BAC->AddElement(fNistManager->FindOrBuildElement(7), 2);

	//Epoxy adhesive with mix ratio 100:25 parts by weight resin/hardener
	// ***don't know density at cold temperature***:
	G4Material* Epoxy=new G4Material("ATEpoxy", density=1.2*g/cm3, ncomponents=2);
	Epoxy->AddMaterial(EpoxyResin, fractionmass=0.8);
	Epoxy->AddMaterial(Epoxy13BAC, fractionmass=0.2);
	//-----------------------------------------------------------------------------

	//just silicon for now
	G4Material* SiPMT = new G4Material("ATSiPMT", density = 2.329*g/cm3, ncomponents = 1);
	SiPMT->AddElement(fNistManager->FindOrBuildElement(14),100.*perCent);    //Si
}
