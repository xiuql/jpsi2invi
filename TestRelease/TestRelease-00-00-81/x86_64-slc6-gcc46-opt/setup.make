----------> uses
# use GaudiInterface GaudiInterface-* External
#   use GaudiKernel *  (no_version_directory)
#     use GaudiPolicy *  (no_version_directory)
#       use LCG_Settings *  (no_version_directory)
#       use Python * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#       use tcmalloc * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.7p3)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use libunwind v* LCG_Interfaces (no_version_directory) (native_version=5c2cade)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#     use Reflex * LCG_Interfaces (no_version_directory)
#       use LCG_Configuration v*  (no_version_directory)
#         use LCG_Platforms *  (no_version_directory)
#       use LCG_Settings v*  (no_version_directory)
#       use ROOT v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#     use Boost * LCG_Interfaces (no_version_directory) (native_version=1.50.0_python2.7)
#       use LCG_Configuration v*  (no_version_directory)
#       use LCG_Settings v*  (no_version_directory)
#       use Python v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#   use GaudiSvc *  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use CLHEP * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.0.4.5)
#     use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#     use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
# use BesCLHEP * External
#   use CLHEP v* LCG_Interfaces (no_version_directory) (native_version=2.0.4.5)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use HepPDT * LCG_Interfaces (no_version_directory) (native_version=2.06.01)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use BesExternalArea BesExternalArea-* External
# use PartPropSvc *  (no_version_directory)
#   use GaudiPolicy *  (no_version_directory)
#   use GaudiKernel *  (no_version_directory)
#   use HepPDT * LCG_Interfaces (no_version_directory) (native_version=2.06.01)
# use CASTOR * LCG_Interfaces (no_version_directory) (native_version=2.1.13-6)
#   use LCG_Configuration v*  (no_version_directory)
#   use LCG_Settings v*  (no_version_directory)
# use ExHelloWorld * BesExamples
#   use BesPolicy BesPolicy-01-* 
#     use BesCxxPolicy BesCxxPolicy-* 
#       use GaudiPolicy v*  (no_version_directory)
#     use BesFortranPolicy BesFortranPolicy-* 
#       use LCG_Settings v*  (no_version_directory)
#   use GaudiInterface GaudiInterface-01-* External
# use RhopiAlg RhopiAlg-00-* Analysis/Physics
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-* External
#   use DstEvent DstEvent-* Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use BesCLHEP BesCLHEP-* External
#     use EventModel EventModel-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use BesCLHEP BesCLHEP-* External
#     use EventModel EventModel-* Event
#     use EvTimeEvent EvTimeEvent-* Event
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use BesCLHEP BesCLHEP-* External
#       use MdcGeomSvc MdcGeomSvc-* Mdc
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use calibUtil * Calibration
#           use GaudiInterface GaudiInterface-01-* External
#           use facilities * Calibration
#             use BesPolicy BesPolicy-* 
#           use xmlBase * Calibration
#             use BesPolicy * 
#             use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#               use LCG_Configuration v*  (no_version_directory)
#               use LCG_Settings v*  (no_version_directory)
#             use facilities * Calibration
#           use rdbModel * Calibration
#             use BesPolicy * 
#             use facilities * Calibration
#             use xmlBase * Calibration
#             use MYSQL * External
#               use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#                 use LCG_Configuration v*  (no_version_directory)
#                 use LCG_Settings v*  (no_version_directory)
#           use BesROOT BesROOT-00-* External
#             use CASTOR v* LCG_Interfaces (no_version_directory) (native_version=2.1.13-6)
#             use ROOT v* LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#               use LCG_Configuration v*  (no_version_directory)
#               use LCG_Settings v*  (no_version_directory)
#               use GCCXML v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=0.9.0_20120309p2)
#                 use LCG_Configuration v*  (no_version_directory)
#                 use LCG_Settings v*  (no_version_directory)
#               use Python v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#               use xrootd v* LCG_Interfaces (no_version_directory) (native_version=3.2.7)
#                 use LCG_Configuration v*  (no_version_directory)
#                 use LCG_Settings v*  (no_version_directory)
#           use DatabaseSvc DatabaseSvc-* Database
#             use BesPolicy BesPolicy-* 
#             use GaudiInterface GaudiInterface-* External
#             use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#             use sqlite * LCG_Interfaces (no_version_directory) (native_version=3070900)
#               use LCG_Configuration v*  (no_version_directory)
#               use LCG_Settings v*  (no_version_directory)
#             use BesROOT * External
#         use CalibData * Calibration
#           use facilities facilities-* Calibration
#           use GaudiInterface * External
#           use BesROOT BesROOT-00-* External
#         use EventModel EventModel-* Event
#       use RelTable RelTable-* Event
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#       use EventModel EventModel-* Event
#       use Identifier Identifier-* DetectorDescription
#         use BesPolicy BesPolicy-* 
#     use MdcRecEvent MdcRecEvent-* Mdc
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use MdcGeomSvc MdcGeomSvc-* Mdc
#       use RelTable RelTable-* Event
#       use EventModel EventModel-* Event
#       use Identifier Identifier-* DetectorDescription
#       use DstEvent DstEvent-* Event
#     use TofRecEvent TofRecEvent-* Tof
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use Identifier Identifier-* DetectorDescription
#       use EventModel EventModel-* Event
#       use DstEvent * Event
#     use EmcRecEventModel EmcRecEventModel-* Emc
#       use BesPolicy BesPolicy-* 
#       use Identifier Identifier-* DetectorDescription
#       use BesCLHEP BesCLHEP-* External
#       use EventModel EventModel-* Event
#       use DstEvent DstEvent-* Event
#       use EmcRecGeoSvc EmcRecGeoSvc-* Emc
#         use BesPolicy BesPolicy-* 
#         use Identifier Identifier-* DetectorDescription
#         use ROOTGeo ROOTGeo-* DetectorDescription
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-* External
#           use BesCLHEP BesCLHEP-* External
#           use BesROOT BesROOT-* External
#           use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#           use GdmlToRoot GdmlToRoot-* External
#             use BesExternalArea BesExternalArea-* External
#             use BesROOT BesROOT-* External
#           use GdmlManagement GdmlManagement-* DetectorDescription
#             use BesExternalArea BesExternalArea-* External
#         use BesCLHEP BesCLHEP-* External
#         use GaudiInterface GaudiInterface-* External
#     use MucRecEvent MucRecEvent-* Muc
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use Identifier Identifier-* DetectorDescription
#       use EventModel EventModel-* Event
#       use ExtEvent ExtEvent-* Event
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use BesCLHEP BesCLHEP-* External
#         use EventModel EventModel-* Event
#         use DstEvent DstEvent-* Event
#       use MucGeomSvc MucGeomSvc-* Muc
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use Identifier Identifier-* DetectorDescription
#         use ROOTGeo ROOTGeo-* DetectorDescription
#         use BesCLHEP BesCLHEP-* External
#         use BesROOT BesROOT-* External
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use GdmlToRoot GdmlToRoot-* External
#         use G4Geo G4Geo-* DetectorDescription
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-* External
#           use BesCLHEP BesCLHEP-* External
#           use BesGeant4 BesGeant4-* External
#             use BesExternalArea BesExternalArea-00-* External
#             use BesCLHEP BesCLHEP-00-* External
#           use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#           use GdmlToG4 GdmlToG4-* External
#             use BesExternalArea BesExternalArea-* External
#             use BesGeant4 BesGeant4-* External
#             use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#           use GdmlManagement GdmlManagement-* DetectorDescription
#           use Identifier Identifier-* DetectorDescription
#           use SimUtil SimUtil-* Simulation/BOOST
#             use BesPolicy BesPolicy-01-* 
#             use BesGeant4 BesGeant4-00-* External
#       use DstEvent DstEvent-* Event
#     use ExtEvent ExtEvent-* Event
#     use DstEvent DstEvent-* Event
#   use McTruth McTruth-* Event
#     use BesPolicy BesPolicy-01-* 
#     use EventModel EventModel-* Event
#     use GaudiInterface GaudiInterface-01-* External
#     use Identifier Identifier-* DetectorDescription
#     use RelTable RelTable-* Event
#   use VertexFit VertexFit-* Analysis
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use BesCLHEP BesCLHEP-* External
#     use MYSQL MYSQL-* External
#     use DstEvent DstEvent-* Event
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use MagneticField MagneticField-* 
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-* External
#       use BesCLHEP * External
#       use BesROOT * External
#       use EventModel EventModel-* Event
#       use rdbModel * Calibration
#       use DatabaseSvc * Database
#       use BesTimerSvc BesTimerSvc-* Utilities
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-* External
#     use EventModel EventModel-* Event
#     use EvtRecEvent EvtRecEvent-* Event
#     use DatabaseSvc DatabaseSvc-* Database
#   use ParticleID ParticleID-* Analysis
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use DstEvent DstEvent-* Event
#     use BesROOT BesROOT-* External
#     use EvtRecEvent EvtRecEvent-* Event
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use TofRecEvent TofRecEvent-* Tof
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use MucRecEvent MucRecEvent-* Muc
#     use ExtEvent ExtEvent-* Event
#   use BesROOT BesROOT-00-* External
# use PipiJpsiAlg PipiJpsiAlg-* Analysis/Physics/PsiPrime
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use TrigEvent TrigEvent-* Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use McTruth McTruth-* Event
#   use MdcRecEvent MdcRecEvent-* Mdc
# use AbsCor AbsCor-* Analysis/PhotonCor
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawEvent RawEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use Identifier Identifier-* DetectorDescription
#       use EventModel EventModel-* Event
#     use EventModel EventModel-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use Identifier Identifier-* DetectorDescription
#   use EmcRecEventModel EmcRecEventModel-* Emc
# use CalibSvc CalibSvc-* Calibration
#   use CalibMySQLCnv * Calibration/CalibSvc
#     use BesPolicy * 
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use GaudiInterface * External
#     use MYSQL MYSQL-00-* External
#     use CalibDataSvc * Calibration/CalibSvc
#       use BesPolicy * 
#       use BesROOT * External
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use DstEvent DstEvent-* Event
#       use EventModel EventModel-* Event
#       use GaudiKernel *  (no_version_directory)
#   use CalibROOTCnv * Calibration/CalibSvc
#     use BesPolicy * 
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use GaudiInterface * External
#     use BesROOT BesROOT-00-* External
#     use EventModel EventModel-* Event
#     use CalibDataSvc * Calibration/CalibSvc
#     use CalibMySQLCnv * Calibration/CalibSvc
#   use CalibDataSvc * Calibration/CalibSvc
#   use CalibXmlCnvSvc * Calibration/CalibSvc
#     use BesPolicy * 
#     use calibUtil * Calibration
#     use GaudiInterface * External
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use CalibDataSvc * Calibration/CalibSvc
#   use CalibTreeCnv * Calibration/CalibSvc
#     use BesPolicy * 
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use MYSQL MYSQL-00-* External
#     use GaudiInterface * External
#     use BesROOT BesROOT-00-* External
#     use DstEvent DstEvent-* Event
#     use EventModel EventModel-* Event
#     use CalibDataSvc * Calibration/CalibSvc
#     use CalibMySQLCnv * Calibration/CalibSvc
#     use DatabaseSvc DatabaseSvc-* Database
# use tofcalgsec tofcalgsec-* Tof
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-* External
#   use BesROOT BesROOT-00-* External
#   use BesCLHEP BesCLHEP-* External
#   use EventModel EventModel-* Event
#   use TofRecEvent TofRecEvent-* Tof
#   use TofRec TofRec-* Reconstruction
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use Identifier Identifier-* DetectorDescription
#     use TofRawEvent TofRawEvent-* Tof
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use RawEvent RawEvent-* Event
#       use EventModel EventModel-* Event
#     use TofRecEvent TofRecEvent-* Tof
#     use TofGeomSvc TofGeomSvc-* Tof
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-* External
#     use TofCaliSvc TofCaliSvc-* Tof
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-* External
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use CalibSvc * Calibration
#       use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#       use MYSQL MYSQL-00-00-* External
#       use BesROOT BesROOT-* External
#     use EvTimeEvent EvTimeEvent-* Event
#     use AsciiDmp AsciiDmp-* Event
#       use BesPolicy BesPolicy-* 
#     use EventModel EventModel-* Event
#     use RawEvent RawEvent-* Event
#     use ExtEvent ExtEvent-* Event
#     use McTruth McTruth-* Event
#     use ReconEvent ReconEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use BesCLHEP BesCLHEP-* External
#       use EventModel EventModel-* Event
#       use ExtEvent ExtEvent-* Event
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#       use BesPolicy BesPolicy-01-* 
#       use RootPolicy RootPolicy-* 
#         use BesPolicy BesPolicy-* 
#         use BesROOT BesROOT-00-* External
#       use GaudiInterface GaudiInterface-01-* External
#       use BesCLHEP BesCLHEP-* External
#       use MdcRawEvent MdcRawEvent-* Mdc
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use TofRawEvent TofRawEvent-* Tof
#       use TofCaliSvc TofCaliSvc-* Tof
#       use TofQCorrSvc TofQCorrSvc-* Tof
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use CalibData * Calibration
#         use EventModel EventModel-* Event
#         use DatabaseSvc DatabaseSvc-* Database
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use MYSQL MYSQL-00-00-* External
#         use BesROOT BesROOT-* External
#       use TofQElecSvc TofQElecSvc-* Tof
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use calibUtil * Calibration
#         use CalibData * Calibration
#         use CalibSvc * Calibration
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use MYSQL MYSQL-00-00-* External
#         use BesROOT BesROOT-* External
#       use EmcRawEvent EmcRawEvent-* Emc
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use RawEvent RawEvent-* Event
#         use EventModel EventModel-* Event
#         use Identifier Identifier-* DetectorDescription
#         use EmcWaveform EmcWaveform-* Emc
#           use BesPolicy BesPolicy-* 
#           use GaudiInterface GaudiInterface-* External
#       use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-* External
#         use CalibData CalibData-* Calibration
#         use CalibDataSvc CalibDataSvc-* Calibration/CalibSvc
#         use CalibROOTCnv CalibROOTCnv-* Calibration/CalibSvc
#         use EmcRecGeoSvc EmcRecGeoSvc-* Emc
#         use EmcGeneralClass EmcGeneralClass-* Emc
#           use BesPolicy BesPolicy-* 
#           use Identifier Identifier-* DetectorDescription
#         use BesCLHEP BesCLHEP-* External
#       use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use CalibData CalibData-* Calibration
#         use CalibSvc CalibSvc-* Calibration
#         use MdcGeomSvc MdcGeomSvc-* Mdc
#         use BesCLHEP BesCLHEP-* External
#       use EvTimeEvent EvTimeEvent-* Event
#       use MrpcCalibSvc MrpcCalibSvc-* Mrpc
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use calibUtil * Calibration
#         use CalibData * Calibration
#         use CalibSvc * Calibration
#         use RawEvent RawEvent-* Event
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use TrigEvent TrigEvent-* Event
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use RootHistCnv v*  (no_version_directory)
#       use GaudiKernel *  (no_version_directory)
#       use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#       use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#   use ExtEvent ExtEvent-* Event
#   use BesGeant4 BesGeant4-00-* External
# use OfflineEventLoopMgr * Control
#   use BesPolicy BesPolicy-* 
#   use BesBoost BesBoost-* External
#     use Boost v* LCG_Interfaces (no_version_directory) (native_version=1.50.0_python2.7)
#   use GaudiInterface GaudiInterface-01-* External
#   use GaudiCoreSvc v*  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#   use GaudiCommonSvc v*  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#     use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#     use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#   use EventModel EventModel-* Event
# use Reconstruction Reconstruction-* 
#   use MdcDedxAlg MdcDedxAlg-* Reconstruction
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use EventModel EventModel-* Event
#     use EvtRecEvent EvtRecEvent-* Event
#     use BesCLHEP BesCLHEP-* External
#     use ReconEvent ReconEvent-* Event
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use DedxCorrecSvc DedxCorrecSvc-* Mdc
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use BesCLHEP BesCLHEP-* External
#       use MdcGeomSvc MdcGeomSvc-* Mdc
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use MagneticField MagneticField-* 
#       use Identifier * DetectorDescription
#       use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#       use MYSQL MYSQL-00-00-* External
#     use EsTimeAlg EsTimeAlg-* Reconstruction
#       use BesPolicy BesPolicy-01-* 
#       use EmcRec EmcRec-* Reconstruction
#         use BesPolicy BesPolicy-* 
#         use Identifier Identifier-* DetectorDescription
#         use GaudiInterface GaudiInterface-01-* External
#         use BesROOT BesROOT-00-* External
#         use EmcRecEventModel EmcRecEventModel-01-* Emc
#         use EmcRecGeoSvc EmcRecGeoSvc-01-* Emc
#         use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#         use RawDataProviderSvc RawDataProviderSvc-* Event
#         use MdcRawEvent MdcRawEvent-* Mdc
#         use TofRawEvent TofRawEvent-* Tof
#         use EmcRawEvent EmcRawEvent-* Emc
#         use EventModel EventModel-* Event
#         use ReconEvent ReconEvent-* Event
#         use McTruth McTruth-* Event
#         use RootHistCnv v*  (no_version_directory)
#       use MdcRawEvent MdcRawEvent-* Mdc
#       use TofRawEvent TofRawEvent-* Tof
#       use EmcRawEvent EmcRawEvent-* Emc
#       use EmcRecEventModel EmcRecEventModel-* Emc
#       use CERNLIB CERNLIB-* External
#         use cernlib v* LCG_Interfaces (no_version_directory) (native_version=2006a)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#           use blas v* LCG_Interfaces (no_version_directory) (native_version=20110419)
#             use LCG_Configuration v*  (no_version_directory)
#             use LCG_Settings v*  (no_version_directory)
#           use lapack v* LCG_Interfaces (no_version_directory) (native_version=3.4.0)
#             use LCG_Configuration v*  (no_version_directory)
#             use LCG_Settings v*  (no_version_directory)
#             use blas v* LCG_Interfaces (no_version_directory) (native_version=20110419)
#         use CASTOR v* LCG_Interfaces (no_version_directory) (native_version=2.1.13-6)
#         use Xt * External (native_version=X11R6)
#       use BesCLHEP BesCLHEP-* External
#       use MdcTables MdcTables-* Mdc
#         use MdcRecEvent MdcRecEvent-* Mdc
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use MdcGeomSvc MdcGeomSvc-* Mdc
#         use Identifier Identifier-* DetectorDescription
#       use MdcGeomSvc MdcGeomSvc-* Mdc
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use EstTofCaliSvc EstTofCaliSvc-* Tof
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use calibUtil * Calibration
#         use CalibData * Calibration
#         use CalibSvc * Calibration
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use MYSQL MYSQL-00-00-* External
#         use BesROOT BesROOT-* External
#       use TofQElecSvc TofQElecSvc-* Tof
#       use TofSim TofSim-* Simulation/BOOST
#         use BesPolicy BesPolicy-01-* 
#         use BesGeant4 BesGeant4-00-* External
#         use GdmlToG4 GdmlToG4-* External
#         use SimUtil SimUtil-* Simulation/BOOST
#         use TruSim TruSim-* Simulation/BOOST
#           use BesPolicy BesPolicy-01-* 
#           use BesGeant4 BesGeant4-* External
#           use BesCLHEP BesCLHEP-* External
#           use GeneratorObject GeneratorObject-* Generator
#             use BesPolicy BesPolicy-* 
#             use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#             use EventModel EventModel-* Event
#         use GaudiInterface GaudiInterface-* External
#         use G4Svc G4Svc-* Simulation
#           use BesCLHEP BesCLHEP-* External
#           use BesPolicy BesPolicy-* 
#           use GaudiInterface GaudiInterface-01-* External
#           use GeneratorObject GeneratorObject-* Generator
#           use BesGeant4 BesGeant4-00-* External
#           use RealizationSvc RealizationSvc-* Simulation/Realization
#             use BesPolicy BesPolicy-01-* 
#             use GaudiInterface GaudiInterface-01-* External
#             use BesCLHEP * External
#             use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#             use EventModel EventModel-* Event
#         use GdmlManagement GdmlManagement-* DetectorDescription
#         use G4Geo G4Geo-* DetectorDescription
#         use TofCaliSvc TofCaliSvc-* Tof
#         use TofSimSvc TofSimSvc-* Tof
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-* External
#           use calibUtil * Calibration
#           use CalibData * Calibration
#           use CalibSvc * Calibration
#           use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#           use MYSQL MYSQL-00-00-* External
#         use TofQElecSvc TofQElecSvc-* Tof
#       use EventModel EventModel-* Event
#       use McTruth McTruth-* Event
#       use RawEvent RawEvent-* Event
#       use RawDataProviderSvc RawDataProviderSvc-* Event
#       use TrigEvent TrigEvent-* Event
#       use EvTimeEvent EvTimeEvent-* Event
#       use MdcFastTrkAlg MdcFastTrkAlg-* Reconstruction
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use ReconEvent ReconEvent-* Event
#         use EventModel EventModel-* Event
#         use MdcGeomSvc MdcGeomSvc-* Mdc
#         use McTruth McTruth-* Event
#         use MdcRecEvent MdcRecEvent-* Mdc
#         use BesCLHEP BesCLHEP-* External
#         use CERNLIB CERNLIB-* External
#         use PartPropSvc *  (no_version_directory)
#         use Identifier Identifier-* DetectorDescription
#         use BesAIDA BesAIDA-* External
#           use AIDA v* LCG_Interfaces (no_version_directory) (native_version=3.2.1)
#             use LCG_Configuration v*  (no_version_directory)
#             use LCG_Settings v*  (no_version_directory)
#         use EvTimeEvent EvTimeEvent-* Event
#         use RawDataProviderSvc RawDataProviderSvc-* Event
#         use TrigEvent TrigEvent-* Event
#         use MagneticField MagneticField-* 
#         use TrackUtil TrackUtil-* Reconstruction
#           use MagneticField MagneticField-* 
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-01-* External
#           use BesCLHEP BesCLHEP-* External
#           use MdcRawEvent MdcRawEvent-* Mdc
#           use MdcRecEvent MdcRecEvent-* Mdc
#           use CERNLIB CERNLIB-* External
#         use RootHistCnv v*  (no_version_directory)
#         use BesTimerSvc BesTimerSvc-00-* Utilities
#       use MdcDedxAlg MdcDedxAlg-* Reconstruction
#       use TrackUtil TrackUtil-* Reconstruction
#       use Identifier Identifier-* DetectorDescription
#       use BesAIDA BesAIDA-* External
#       use BesTimerSvc BesTimerSvc-00-* Utilities
#       use PartPropSvc *  (no_version_directory)
#       use RootHistCnv v*  (no_version_directory)
#       use BesROOT BesROOT-* External
#       use RawEvent * Event
#       use MagneticField MagneticField-* 
#       use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#       use MdcUtilitySvc MdcUtilitySvc-* Mdc/MdcCheckUtil
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use MagneticField MagneticField-* 
#         use MdcRecEvent MdcRecEvent-* Mdc
#         use MdcGeomSvc MdcGeomSvc-* Mdc
#         use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-* External
#           use Identifier Identifier-* DetectorDescription
#           use MdcGeomSvc MdcGeomSvc-* Mdc
#           use BesCLHEP BesCLHEP-* External
#         use TrkBase TrkBase-* Reconstruction/MdcPatRec
#           use BesPolicy BesPolicy-01-* 
#           use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#           use MdcRecoUtil MdcRecoUtil-* Reconstruction/MdcPatRec
#             use BesPolicy BesPolicy-* 
#             use BesCLHEP BesCLHEP-* External
#           use ProxyDict ProxyDict-* Reconstruction/MdcPatRec
#             use BesPolicy BesPolicy-01-* 
#           use ProbTools ProbTools-* Reconstruction/MdcPatRec
#             use BesPolicy BesPolicy-01-* 
#           use BField BField-* Reconstruction/MdcPatRec
#             use BesPolicy BesPolicy-* 
#             use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#             use MagneticField MagneticField-* 
#           use BesCLHEP BesCLHEP-* External
#           use CERNLIB CERNLIB-* External
#           use BesBoost BesBoost-* External
#           use MdcRecEvent MdcRecEvent-* Mdc
#     use DedxCurSvc DedxCurSvc-* Mdc
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-* External
#       use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#       use MYSQL MYSQL-00-00-* External
#       use EventModel EventModel-* Event
#       use DatabaseSvc DatabaseSvc-* Database
#       use BesROOT BesROOT-* External
#     use CERNLIB CERNLIB-* External
#     use BesAIDA BesAIDA-* External
#     use RootHistCnv v*  (no_version_directory)
#   use MucRecAlg MucRecAlg-* Reconstruction
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use MucRawEvent MucRawEvent-* Muc
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use RawEvent RawEvent-* Event
#       use EventModel EventModel-* Event
#     use EventModel EventModel-* Event
#     use ExtEvent ExtEvent-* Event
#     use ReconEvent ReconEvent-* Event
#     use Identifier Identifier-* DetectorDescription
#     use PartPropSvc *  (no_version_directory)
#     use MucGeomSvc MucGeomSvc-* Muc
#     use MucRecEvent MucRecEvent-* Muc
#     use McTruth McTruth-* Event
#     use EmcRecEventModel EmcRecEventModel-01-* Emc
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use RootHistCnv v*  (no_version_directory)
#   use TofRec TofRec-* Reconstruction
#   use MdcTrkRecon MdcTrkRecon-* Reconstruction/MdcPatRec
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use EventModel EventModel-* Event
#     use McTruth McTruth-* Event
#     use BesCLHEP BesCLHEP-* External
#     use BesBoost BesBoost-* External
#     use BesAIDA BesAIDA-* External
#     use MdcData MdcData-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use EventModel EventModel-* Event
#       use RawEvent RawEvent-* Event
#       use MdcRawEvent MdcRawEvent-* Mdc
#       use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#       use Identifier Identifier-* DetectorDescription
#       use TrkBase TrkBase-* Reconstruction/MdcPatRec
#       use MdcRecoUtil MdcRecoUtil-* Reconstruction/MdcPatRec
#     use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#     use BField BField-* Reconstruction/MdcPatRec
#     use TrkFitter TrkFitter-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-01-* 
#       use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#       use BField BField-* Reconstruction/MdcPatRec
#       use BesCLHEP BesCLHEP-* External
#       use TrkBase TrkBase-* Reconstruction/MdcPatRec
#     use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#     use MagneticField MagneticField-* 
#     use RootHistCnv v*  (no_version_directory)
#     use EvTimeEvent EvTimeEvent-* Event
#     use ReconEvent ReconEvent-* Event
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#     use MdcPrintSvc MdcPrintSvc-* Mdc/MdcCheckUtil
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-* External
#       use EventModel EventModel-* Event
#       use McTruth McTruth-* Event
#       use RawDataProviderSvc RawDataProviderSvc-* Event
#       use EvtRecEvent EvtRecEvent-* Event
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use MdcRawEvent MdcRawEvent-* Mdc
#   use EmcRec EmcRec-* Reconstruction
#   use MdcFastTrkAlg MdcFastTrkAlg-* Reconstruction
#   use TrkExtAlg TrkExtAlg-* Reconstruction
#     use BesPolicy BesPolicy-01-* 
#     use BesCLHEP BesCLHEP-* External
#     use GdmlToG4 GdmlToG4-* External
#     use BesGeant4 BesGeant4-00-* External
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcRecEvent MdcRecEvent* Mdc
#     use EventModel EventModel-* Event
#     use ExtEvent ExtEvent-* Event
#     use TofSim TofSim-* Simulation/BOOST
#     use SimUtil SimUtil-* Simulation/BOOST
#     use G4Geo G4Geo-* DetectorDescription
#     use TrackUtil TrackUtil-* Reconstruction
#     use MagneticField MagneticField-* 
#     use MucRawEvent MucRawEvent-* Muc
#     use MucGeomSvc MucGeomSvc-* Muc
#     use MucRecEvent MucRecEvent-* Muc
#     use MucCalibAlg MucCalibAlg-* Muc
#       use BesPolicy BesPolicy-* 
#       use Identifier Identifier-* DetectorDescription
#       use GaudiInterface GaudiInterface-* External
#       use McTruth McTruth-* Event
#       use EventModel EventModel-* Event
#       use DstEvent DstEvent-* Event
#       use EvTimeEvent EvTimeEvent-* Event
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use TofRecEvent TofRecEvent-* Tof
#       use EmcRecEventModel EmcRecEventModel-* Emc
#       use MucRawEvent MucRawEvent-* Muc
#       use MucRecEvent MucRecEvent-* Muc
#       use RootHistCnv v*  (no_version_directory)
#   use EsTimeAlg EsTimeAlg-* Reconstruction
# use MdcDedxAlg MdcDedxAlg-* Reconstruction
# use MucRecAlg MucRecAlg-* Reconstruction
# use TofRec TofRec-* Reconstruction
# use TofEnergyRec TofEnergyRec-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use Identifier Identifier-* DetectorDescription
#   use TofRawEvent TofRawEvent-* Tof
#   use TofRecEvent TofRecEvent-* Tof
#   use TofGeomSvc TofGeomSvc-* Tof
#   use TofCaliSvc TofCaliSvc-* Tof
#   use TofQCorrSvc TofQCorrSvc-* Tof
#   use TofQElecSvc TofQElecSvc-* Tof
#   use TofSim TofSim-* Simulation/BOOST
#   use TofRec TofRec-* Reconstruction
#   use MrpcCalibSvc MrpcCalibSvc-* Mrpc
#   use EvTimeEvent EvTimeEvent-* Event
#   use AsciiDmp AsciiDmp-* Event
#   use EventModel EventModel-* Event
#   use RawEvent RawEvent-* Event
#   use ExtEvent ExtEvent-* Event
#   use McTruth McTruth-* Event
#   use ReconEvent ReconEvent-* Event
#   use RawDataProviderSvc RawDataProviderSvc-* Event
#   use RootHistCnv v*  (no_version_directory)
# use MdcTrkRecon MdcTrkRecon-* Reconstruction/MdcPatRec
# use MdcxReco MdcxReco-* Reconstruction/MdcPatRec
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-* External
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use EventModel EventModel-* Event
#   use McTruth McTruth-* Event
#   use BesCLHEP BesCLHEP-* External
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use EvTimeEvent EvTimeEvent-* Event
#   use Identifier Identifier-* DetectorDescription
#   use BesAIDA BesAIDA-* External
#   use MdcData MdcData-* Reconstruction/MdcPatRec
#   use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#   use BField BField-* Reconstruction/MdcPatRec
#   use TrkFitter TrkFitter-* Reconstruction/MdcPatRec
#   use TrkBase TrkBase-* Reconstruction/MdcPatRec
#   use MdcRecoUtil MdcRecoUtil-* Reconstruction/MdcPatRec
#   use MdcTrkRecon MdcTrkRecon-* Reconstruction/MdcPatRec
#   use TrigEvent TrigEvent-* Event
#   use MdcUtilitySvc MdcUtilitySvc-* Mdc/MdcCheckUtil
#   use MdcPrintSvc MdcPrintSvc-* Mdc/MdcCheckUtil
#   use RootHistCnv v*  (no_version_directory)
# use EmcRec EmcRec-* Reconstruction
# use EmcTimeRec EmcTimeRec-* Reconstruction
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use EventModel EventModel-* Event
#   use EvTimeEvent EvTimeEvent-* Event
#   use EmcRecEventModel EmcRecEventModel-01-* Emc
# use MdcFastTrkAlg MdcFastTrkAlg-* Reconstruction
# use TrkExtAlg TrkExtAlg-* Reconstruction
# use EsTimeAlg EsTimeAlg-* Reconstruction
# use KalFitAlg KalFitAlg-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-* External
#   use CERNLIB CERNLIB-* External
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use BesCLHEP BesCLHEP-* External
#   use EventModel EventModel-* Event
#   use EvTimeEvent EvTimeEvent-* Event
#   use MdcTables MdcTables-* Mdc
#   use ReconEvent ReconEvent-* Event
#   use McTruth McTruth-* Event
#   use G4Svc G4Svc-00-* Simulation
#   use G4Geo * DetectorDescription
#   use RootHistCnv v*  (no_version_directory)
#   use BesROOT BesROOT-* External
#   use BesAIDA * External
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use MagneticField MagneticField-* 
#   use PartPropSvc *  (no_version_directory)
#   use VertexFit VertexFit-* Analysis
# use TrkReco TrkReco-* Reconstruction
#   use RawDataProviderSvc RawDataProviderSvc-* Event
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use RawEvent RawEvent-* Event
#   use EvTimeEvent EvTimeEvent-* Event
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use CERNLIB CERNLIB-* External
#   use BesCLHEP BesCLHEP-* External
#   use EventModel EventModel-* Event
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use MdcTables MdcTables-* Mdc
#   use McTruth McTruth-* Event
#   use Identifier Identifier-* DetectorDescription
#   use RawDataProviderSvc RawDataProviderSvc-* Event
#   use TrackUtil TrackUtil-* Reconstruction
#   use RootHistCnv v*  (no_version_directory)
#   use BesAIDA * External
#   use BesTimerSvc BesTimerSvc-00-* Utilities
#   use G4Svc G4Svc-00-* Simulation
#   use G4Geo * DetectorDescription
#   use ReconEvent ReconEvent-* Event
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
# use EventAssembly EventAssembly-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use EvtRecEvent EvtRecEvent-* Event
#   use EventModel EventModel-* Event
#   use EvTimeEvent EvTimeEvent-* Event
# use PrimaryVertexAlg PrimaryVertexAlg-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use McTruth McTruth-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use EvtRecEvent EvtRecEvent-* Event
#   use ReconEvent ReconEvent-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
# use VeeVertexAlg VeeVertexAlg-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use McTruth McTruth-* Event
#   use ReconEvent ReconEvent-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-* External
# use BeamParamsAlg BeamParamsAlg-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use McTruth McTruth-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use EvtRecEvent EvtRecEvent-* Event
#   use ReconEvent ReconEvent-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
# use DTagSkim DTagSkim-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use BesAIDA BesAIDA-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#   use EmcRawEvent EmcRawEvent-* Emc
#   use DTagAlg DTagAlg-* Reconstruction
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use DstEvent DstEvent-* Event
#     use EventModel EventModel-* Event
#     use EvtRecEvent EvtRecEvent-* Event
#     use DecayChain DecayChain-* Event
#       use BesPolicy BesPolicy-* 
#     use BesDChain BesDChain-* Event
#       use BesPolicy BesPolicy-* 
#       use BesCLHEP BesCLHEP-* External
#       use EvtRecEvent EvtRecEvent-* Event
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use EmcRecEventModel EmcRecEventModel-* Emc
#       use VertexFit VertexFit-* Analysis
#       use DecayChain DecayChain-* Event
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use McTruth McTruth-* Event
#     use ReconEvent ReconEvent-* Event
#     use ParticleID ParticleID-* Analysis
#     use VertexFit VertexFit-* Analysis
#     use BesROOT BesROOT-* External
#     use SimplePIDSvc SimplePIDSvc-* Utilities
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use EvtRecEvent EvtRecEvent-* Event
#     use DatabaseSvc DatabaseSvc-* Database
#   use RootIO RootIO-* Event
#     use GaudiKernel *  (no_version_directory)
#     use Gaudi *  (no_version_directory)
#       use GaudiSys *  (no_version_directory)
#         use GaudiKernel *  (no_version_directory)
#         use GaudiCoreSvc *  (no_auto_imports) (no_version_directory)
#         use GaudiCommonSvc *  (no_auto_imports) (no_version_directory)
#         use GaudiAud *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#         use GaudiAlg *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#           use GaudiUtils *  (no_version_directory)
#             use GaudiKernel *  (no_version_directory)
#             use ROOT * LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#             use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#             use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#             use uuid * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.42)
#               use LCG_Configuration v*  (no_version_directory)
#               use LCG_Settings v*  (no_version_directory)
#             use Reflex * LCG_Interfaces (no_auto_imports) (no_version_directory)
#             use XercesC * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.1.1p1)
#           use ROOT * LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#           use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#           use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#         use GaudiPython *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#           use GaudiAlg *  (no_version_directory)
#           use GaudiUtils *  (no_version_directory)
#           use AIDA * LCG_Interfaces (no_version_directory) (native_version=3.2.1)
#           use CLHEP * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.0.4.5)
#           use Python * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#           use Reflex * LCG_Interfaces (no_auto_imports) (no_version_directory)
#       use pytools * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.8_python2.7)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use Python v* LCG_Interfaces (no_version_directory) (native_version=2.7.3)
#         use mysql v* LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#     use ReconEvent ReconEvent-* Event
#     use RootCnvSvc RootCnvSvc-* Event
#       use BesCLHEP BesCLHEP-* External
#       use BesPolicy BesPolicy-01-* 
#       use DataInfoSvc DataInfoSvc-* Control
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#       use BesROOT BesROOT-00-* External
#       use GaudiInterface GaudiInterface-* External
#       use Identifier Identifier-* DetectorDescription
#       use RawEvent RawEvent-* Event
#       use RootEventData RootEventData-* Event
#         use RootPolicy RootPolicy-* 
#         use BesROOT BesROOT-* External
#         use MucRecEvent MucRecEvent-* Muc
#         use Identifier Identifier-* DetectorDescription
#       use MdcRawEvent MdcRawEvent-* Mdc
#       use EventModel EventModel-* Event
#       use EmcRawEvent EmcRawEvent-* Emc
#       use TofRawEvent TofRawEvent-* Tof
#       use TofRecEvent TofRecEvent-* Tof
#       use MucRawEvent MucRawEvent-* Muc
#       use MucRecEvent MucRecEvent-* Muc
#       use EmcRecEventModel EmcRecEventModel-* Emc
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use EvTimeEvent EvTimeEvent-* Event
#       use ZddEvent ZddEvent-* Event
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-* External
#       use EvtRecEvent EvtRecEvent-* Event
#       use EventNavigator EventNavigator-* Event
#         use BesPolicy * 
#         use GaudiInterface * External
#         use McTruth McTruth-* Event
#         use EmcRecEventModel * Emc
#         use MdcRecEvent * Mdc
#         use MdcRawEvent * Mdc
#         use MucRecEvent * Muc
#         use TofRecEvent * Tof
#       use McTruth McTruth-* Event
#       use ExtEvent ExtEvent-* Event
#       use DstEvent DstEvent-* Event
#       use ReconEvent ReconEvent-* Event
#       use TrigEvent TrigEvent-* Event
#       use HltEvent HltEvent-* Event
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use RawEvent RawEvent-* Event
#         use EventModel EventModel-* Event
#         use Identifier Identifier-* DetectorDescription
#       use LumiDigi LumiDigi-* Event
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use RawEvent RawEvent-* Event
#         use TofRawEvent TofRawEvent-* Tof
#     use BesROOT BesROOT-00-* External
#     use BesCLHEP BesCLHEP-* External
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use DstEvent DstEvent-* Event
#     use EventModel EventModel-* Event
#     use McTruth McTruth-* Event
#     use Identifier Identifier-* DetectorDescription
#     use RootPolicy RootPolicy-* 
#     use RootEventData RootEventData-* Event
# use LumTauAlg LumTauAlg-* 
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use McTruth McTruth-* Event
#   use RootHistCnv v*  (no_version_directory)
# use T0Dummy T0Dummy-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use EventModel EventModel-* Event
#   use McTruth McTruth-* Event
#   use EvTimeEvent EvTimeEvent-* Event
# use DiGamAlg DiGamAlg-* DQA
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use TofRecEvent TofRecEvent-* Tof
#   use EmcRecEventModel EmcRecEventModel-* Emc
#   use MucRecEvent MucRecEvent-* Muc
#   use ExtEvent ExtEvent-* Event
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#   use McTruth McTruth-* Event
#   use VertexFit VertexFit-* Analysis
#   use ParticleID ParticleID-* Analysis
#   use RootHistCnv v*  (no_version_directory)
#   use BesROOT BesROOT-00-* External
# use BesSim BesSim-* Simulation/BOOST
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use BesGeant4 BesGeant4-00-* External
#   use GdmlToG4 GdmlToG4-* External
#   use G4Svc G4Svc-00-* Simulation
#   use BesServices BesServices-* Control
#     use BesPolicy BesPolicy-* 
#     use BesKernel BesKernel-* Control
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#     use BesBoost BesBoost-* External
#     use CLHEP * LCG_Interfaces (no_version_directory) (native_version=2.0.4.5)
#     use GaudiInterface GaudiInterface-* External
#     use GaudiSvc v*  (no_version_directory)
#   use BesROOT BesROOT-* External
#   use RootPolicy RootPolicy-* 
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use TofRawEvent TofRawEvent-* Tof
#   use EmcRawEvent EmcRawEvent-* Emc
#   use MucRawEvent MucRawEvent-* Muc
#   use Identifier Identifier-* DetectorDescription
#   use RawEvent RawEvent-* Event
#   use McTruth McTruth-* Event
#   use AsciiDmp AsciiDmp-* Event
#   use MagneticField MagneticField-* 
#   use RealizationSvc RealizationSvc-* Simulation/Realization
#   use RootEventData RootEventData-* Event
#   use SimUtil SimUtil-* Simulation/BOOST
#   use TruSim TruSim-* Simulation/BOOST
#   use MdcSim MdcSim-* Simulation/BOOST
#     use BesPolicy BesPolicy-01-* 
#     use BesGeant4 BesGeant4-00-* External
#     use BesCLHEP BesCLHEP-* External
#     use GdmlToG4 GdmlToG4-* External
#     use SimUtil SimUtil-* Simulation/BOOST
#     use TruSim TruSim-* Simulation/BOOST
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#     use MdcTunningSvc MdcTunningSvc-* Simulation/BOOST
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use EventModel EventModel-* Event
#       use DatabaseSvc * Database
#     use G4Svc G4Svc-00-* Simulation
#     use BesROOT * External
#     use G4Geo G4Geo-* DetectorDescription
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use DedxCurSvc DedxCurSvc-* Mdc
#   use TofSim TofSim-* Simulation/BOOST
#   use EmcSim EmcSim-* Simulation/BOOST
#     use BesPolicy BesPolicy-01-* 
#     use BesGeant4 BesGeant4-00-* External
#     use GdmlToG4 GdmlToG4-* External
#     use SimUtil SimUtil-* Simulation/BOOST
#     use TruSim TruSim-* Simulation/BOOST
#     use GaudiInterface GaudiInterface-* External
#     use G4Svc G4Svc-* Simulation
#     use G4Geo G4Geo-* DetectorDescription
#     use Identifier Identifier-* DetectorDescription
#     use EmcGeneralClass EmcGeneralClass-* Emc
#     use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#   use MucSim MucSim-* Simulation/BOOST
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use BesGeant4 BesGeant4-00-* External
#     use GdmlToG4 GdmlToG4-* External
#     use SimUtil SimUtil-* Simulation/BOOST
#     use TruSim TruSim-* Simulation/BOOST
#     use G4Geo G4Geo-* DetectorDescription
#     use BesROOT BesROOT-* External
#     use G4Svc G4Svc-* Simulation
#     use MucCalibConstSvc MucCalibConstSvc-* Muc
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use CalibSvc CalibSvc-* Calibration
#       use CalibData CalibData-* Calibration
#       use CalibDataSvc CalibDataSvc-* Calibration/CalibSvc
#     use MucCalibAlg MucCalibAlg-* Muc
#   use GenSim GenSim-* Simulation/BOOST
#     use BesPolicy BesPolicy-01-* 
#     use BesGeant4 BesGeant4-00-* External
#     use BesROOT BesROOT-* External
#   use PhySim PhySim-* Simulation/BOOST
#     use BesPolicy BesPolicy-01-* 
#     use BesGeant4 BesGeant4-00-* External
#   use EventNavigator EventNavigator-* Event
# use MdcTunningSvc MdcTunningSvc-* Simulation/BOOST
# use BesEventMixer BesEventMixer-* Simulation
#   use BesPolicy * 
#   use GaudiInterface * External
#   use BesCLHEP * External
#   use DataInfoSvc DataInfoSvc-* Control
#   use IRawFile IRawFile-* Event
#     use BesPolicy BesPolicy-01-* 
#   use RawFile RawFile-* Event
#     use BesPolicy BesPolicy-* 
#     use BesExternalArea BesExternalArea-* External
#     use IRawFile IRawFile-* Event
#     use XmlRpcC XmlRpcC-* External
#       use BesExternalArea BesExternalArea-* External
#   use RawDataCnv * Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawFile RawFile-* Event
#     use EventModel EventModel-* Event
#     use RawDataCnvBase RawDataCnvBase-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use EventModel EventModel-* Event
#     use RawEvent RawEvent-00-* Event
#     use McTruth McTruth-* Event
#     use Identifier Identifier-* DetectorDescription
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use MucRawEvent MucRawEvent-* Muc
#     use TrigEvent TrigEvent-* Event
#     use HltEvent HltEvent-* Event
#     use LumiDigi LumiDigi-* Event
#     use ZddEvent ZddEvent-* Event
#     use eformat eformat-* Event
#       use BesPolicy BesPolicy-* 
#       use ers ers-* Event
#         use BesPolicy BesPolicy-* 
#     use TofSim TofSim-* Simulation/BOOST
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#   use RawEvent * Event
#   use RootCnvSvc * Event
#   use HltEvent HltEvent-* Event
#   use TrigEvent TrigEvent-* Event
#   use RealizationSvc RealizationSvc-* Simulation/Realization
#   use EventModel EventModel-* Event
#   use BesTimerSvc BesTimerSvc-* Utilities
#   use BesServices BesServices-* Control
# use Trigger Trigger-* Trigger
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use BesROOT * External
#   use BesCLHEP * External
#   use BesAIDA BesAIDA-* External
#   use RealizationSvc RealizationSvc-* Simulation/Realization
#   use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#   use RawDataProviderSvc RawDataProviderSvc-* Event
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use TofRawEvent TofRawEvent-* Tof
#   use EmcRawEvent EmcRawEvent-* Emc
#   use MucRawEvent MucRawEvent-* Muc
#   use EmcWaveform EmcWaveform-* Emc
#   use TrigEvent TrigEvent-* Event
#   use McTruth McTruth-* Event
#   use RawEvent RawEvent-* Event
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#   use Identifier Identifier-* DetectorDescription
# use TrigMakerAlg TrigMakerAlg-* Event
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use TrigEvent TrigEvent-* Event
#   use EventModel EventModel-* Event
# use EventFilter EventFilter-* 
#   use HltTools HltTools-* EventFilter/HltUtilities
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#   use HltDataTypes HltDataTypes-* EventFilter/HltUtilities
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#   use HltStore HltStore-* EventFilter/HltUtilities
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#   use HltAlgorithms HltAlgorithms-* EventFilter/HltUtilities
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use HltStore HltStore-* EventFilter/HltUtilities
#     use HltDataTypes HltDataTypes-* EventFilter/HltUtilities
#     use RawEvent RawEvent-* Event
#     use EmcRawEvent EmcRawEvent-* Emc
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use MucRawEvent MucRawEvent-* Muc
#     use HltEvent HltEvent-* Event
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use EventModel EventModel-* Event
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#     use RootHistCnv v*  (no_version_directory)
#   use HltSteerData HltSteerData-* EventFilter/HltProcessor
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use HltTools HltTools-* EventFilter/HltUtilities
#     use HltStore HltStore-* EventFilter/HltUtilities
#     use HltDataTypes HltDataTypes-* EventFilter/HltUtilities
#   use HltConfig HltConfig-* EventFilter/HltProcessor
#     use BesPolicy BesPolicy-* 
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use GaudiInterface GaudiInterface-* External
#     use DatabaseSvc DatabaseSvc-* Database
#     use EventModel EventModel-* Event
#     use HltSteerData HltSteerData-* EventFilter/HltProcessor
#     use HltTools HltTools-* EventFilter/HltUtilities
#   use HltSteering HltSteering-* EventFilter/HltProcessor
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use HltDataTypes HltDataTypes-* EventFilter/HltUtilities
#     use HltStore HltStore-* EventFilter/HltUtilities
#     use HltAlgorithms HltAlgorithms-* EventFilter/HltUtilities
#     use HltSteerData HltSteerData-* EventFilter/HltProcessor
#     use HltConfig HltConfig-* EventFilter/HltProcessor
#   use HltEvent HltEvent-* Event
#   use TrigEvent TrigEvent-* Event
#   use TofRawEvent TofRawEvent-* Tof
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use EmcRawEvent EmcRawEvent-* Emc
#   use MucRawEvent MucRawEvent-* Muc
#   use efhlt efhlt-* EventFilter/OnlineJointer
#   use EFServices EFServices-* EventFilter/OnlineJointer
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use GaudiSvc v*  (no_version_directory)
#     use RawFile RawFile-* Event
#     use EventModel EventModel-* Event
#     use RawDataCnvBase RawDataCnvBase-* Event
#     use RawEvent RawEvent-00-* Event
#     use McTruth McTruth-* Event
#     use Identifier Identifier-* DetectorDescription
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use MucRawEvent MucRawEvent-* Muc
#     use TrigEvent TrigEvent-* Event
#     use HltEvent HltEvent-* Event
#     use LumiDigi LumiDigi-* Event
#     use RawDataCnv RawDataCnv-* Event
#     use HltDataTypes HltDataTypes-* EventFilter/HltUtilities
#     use HltStore HltStore-* EventFilter/HltUtilities
#   use ESController ESController-* EventFilter/OnlineJointer
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use HltTools HltTools-* EventFilter/HltUtilities
#     use EFServices EFServices-* EventFilter/OnlineJointer
#     use efhlt efhlt-* EventFilter/OnlineJointer
# use HltMakerAlg HltMakerAlg-* Event
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use HltEvent HltEvent-* Event
#   use EventModel EventModel-* Event
#   use Identifier Identifier-* DetectorDescription
#   use TofRawEvent TofRawEvent-* Tof
#   use DstMakerAlg DstMakerAlg-* Event
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use TofRecEvent TofRecEvent-* Tof
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use MucRecEvent MucRecEvent-* Muc
#     use ExtEvent ExtEvent-* Event
#     use BesCLHEP BesCLHEP-* External
#     use BesROOT BesROOT-* External
#     use DstEvent DstEvent-* Event
#     use EventModel EventModel-* Event
# use EventPreSelect EventPreSelect-* EvtPreSelect
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#   use EmcRawEvent EmcRawEvent-* Emc
#   use MucRawEvent MucRawEvent-* Muc
#   use MdcRawEvent * Mdc
# use EventWriter EventWriter-* EvtPreSelect
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use GaudiSvc *  (no_version_directory)
#   use RootPolicy RootPolicy-* 
#   use BesROOT BesROOT-* External
#   use BesCLHEP BesCLHEP-* External
#   use EvtRecEvent EvtRecEvent-* Event
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use TofRecEvent TofRecEvent-* Tof
#   use EmcRecEventModel EmcRecEventModel-* Emc
#   use MucRecEvent MucRecEvent-* Muc
#   use EmcRec EmcRec-* Reconstruction
#   use RootCnvSvc RootCnvSvc-* Event
#   use EventModel EventModel-* Event
#   use ExtEvent ExtEvent-* Event
#   use DstEvent DstEvent-* Event
#   use RootEventData RootEventData-* Event
#   use DataInfoSvc DataInfoSvc-* Control
#   use DistBossUtil DistBossUtil-* DistBoss
#     use BesPolicy BesPolicy-* 
#     use IRawFile IRawFile-* Event
#   use NetDataWriter NetDataWriter-* DistBoss
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use DimSetup DimSetup-* DistBoss
#       use BesExternalArea BesExternalArea-* External
#     use DistBossUtil DistBossUtil-* DistBoss
#     use ClientErrHandler ClientErrHandler-* DistBoss
#       use BesPolicy BesPolicy-* 
#       use DimSetup DimSetup-* DistBoss
#       use DistBossUtil DistBossUtil-* DistBoss
# use CalibEventSelect CalibEventSelect-* EvtPreSelect
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use BesAIDA BesAIDA-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use ParticleID ParticleID-* Analysis
#   use DTagAlg DTagAlg-* Reconstruction
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#   use EmcRawEvent EmcRawEvent-* Emc
#   use RootIO RootIO-* Event
# use DimuPreSelect DimuPreSelect-* EvtPreSelect
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use MucRawEvent MucRawEvent-* Muc
#   use ParticleID ParticleID-* Analysis
#   use VertexFit VertexFit-* Analysis
#   use BesROOT BesROOT-00-* External
#   use EmcRecEventModel EmcRecEventModel-01-* Emc
#   use RootHistCnv v*  (no_version_directory)
# use Gam4pikpAlg Gam4pikpAlg-* Analysis/ControlSamples
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-* External
#   use McTruth McTruth-* Event
#   use DstEvent DstEvent-* Event
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#   use HltEvent HltEvent-* Event
#   use VertexFit VertexFit-* Analysis
#   use ParticleID ParticleID-* Analysis
#   use BesROOT BesROOT-* External
# use RootIO RootIO-* Event
# use RawDataCnv RawDataCnv-* Event
# use RawDataProviderSvc RawDataProviderSvc-* Event
# use MagneticField MagneticField-* 
# use BesVisAlg * EventDisplay
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use RootPolicy RootPolicy-* 
#   use BesROOT BesROOT-* External
#   use BesCLHEP BesCLHEP-* External
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use TofRecEvent TofRecEvent-* Tof
#   use EmcRecEventModel EmcRecEventModel-* Emc
#   use MucRecEvent MucRecEvent-* Muc
#   use RootCnvSvc RootCnvSvc-* Event
#   use EventModel EventModel-* Event
#   use ExtEvent ExtEvent-* Event
#   use DstEvent DstEvent-* Event
#   use BesVisLib BesVisLib-* EventDisplay
#     use BesPolicy BesPolicy-01-* 
#     use RootPolicy RootPolicy-* 
#     use BesROOT BesROOT-* External
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use BesCLHEP BesCLHEP-* External
#     use GdmlToRoot GdmlToRoot-* External
#     use RootEventData RootEventData-* Event
#     use Identifier Identifier-* DetectorDescription
#     use RawEvent RawEvent-* Event
#   use BesVisClient BesVisClient-* EventDisplay
#     use BesPolicy BesPolicy-* 
#     use RootPolicy RootPolicy-* 
#     use BesROOT BesROOT-* External
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use BesCLHEP BesCLHEP-* External
#     use GdmlToRoot GdmlToRoot-* External
#     use RootEventData * Event
#     use Identifier Identifier-* DetectorDescription
#     use BesVisLib BesVisLib-* EventDisplay
# use PartPropSvc *  (no_version_directory)
# use KKMC KKMC-* Generator
#   use BesPolicy BesPolicy-* 
#   use BesFortranPolicy BesFortranPolicy-* 
#   use BesServices BesServices-* Control
#   use CERNLIB CERNLIB-* External
#   use GeneratorObject GeneratorObject-* Generator
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use DatabaseSvc DatabaseSvc-* Database
# use BesEvtGen BesEvtGen-* Generator
#   use BesPolicy BesPolicy-* 
#   use BesServices BesServices-* Control
#   use CERNLIB CERNLIB-* External
#   use GeneratorObject GeneratorObject-* Generator
#   use BesEvtGen BesEvtGen-* Generator
#   use PartPropSvc *  (no_version_directory)
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use DataInfoSvc DataInfoSvc-* Control
#   use BesROOT BesROOT-* External
# use Bhwide Bhwide-* Generator
#   use BesPolicy BesPolicy-* 
#   use BesKernel * Control
#   use CERNLIB CERNLIB-* External
#   use GeneratorObject GeneratorObject-* Generator
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use BesCLHEP BesCLHEP-* External
# use Bhagen Bhagen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#     use BesPolicy BesPolicy-* 
#     use PartPropSvc *  (no_version_directory)
#     use GaudiInterface GaudiInterface-* External
#     use GeneratorObject GeneratorObject-* Generator
#     use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#     use BesServices BesServices-* Control
#   use GeneratorUtil GeneratorUtil-00-* Generator
#     use BesPolicy BesPolicy-* 
#   use bhagen bhagen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#       use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use BesGenModule BesGenModule-* Generator
#     use McTruth McTruth-* Event
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Ddgen Ddgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use ddgen ddgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Ddprod Ddprod-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use ddprod ddprod-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Dsdgen Dsdgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use dsdgen dsdgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Dssgen Dssgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use dssgen dssgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Epscat Epscat-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use epscat epscat-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Fff Fff-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use fff fff-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Ffgen Ffgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use ffgen ffgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Fsfgen Fsfgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use fsfgen fsfgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
# use Gamma2 Gamma2-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use gamma2 gamma2-* External/GENBES
#     use BesExternalArea BesExternalArea-* External
#     use CERNLIB CERNLIB-01-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Howl Howl-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use howl howl-* External/GENBES (native_version=genbes-00-00-11)
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Koralbe Koralbe-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use koralbe koralbe-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use KK2f KK2f-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use kk2f kk2f-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Kstark Kstark-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use kstark kstark-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Lund Lund-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use lund lund-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Lund_crm Lund_crm-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use lund_crm lund_crm-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Mugen Mugen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use mugen mugen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use P2bb P2bb-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use p2bb p2bb-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use P2epem P2epem-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use p2epem p2epem-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use P2mumu P2mumu-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use p2mumu p2mumu-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Ppgen Ppgen-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use ppgen ppgen-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Radee Radee-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use radee radee-* External/GENBES (native_version=genbes-00-00-11)
#     use CERNLIB CERNLIB-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Radgg Radgg-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use radgg radgg-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Radmu Radmu-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use radmu radmu-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Rhopi Rhopi-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use rhopi rhopi-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Sagerx Sagerx-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use sagerx sagerx-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Tauprd Tauprd-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use tauprd tauprd-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Tester Tester-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use tester tester-* External/GENBES (native_version=genbes-00-00-11)
#     use CERNLIB CERNLIB-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
# use Twogam Twogam-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use twogam twogam-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use V2llg V2llg-* Generator/BesGenInterface
#   use BesPolicy BesPolicy-* 
#   use GeneratorModule GeneratorModule-* Generator
#   use GeneratorUtil GeneratorUtil-00-* Generator
#   use v2llg v2llg-* External/GENBES
#     use CERNLIB CERNLIB-01-* External
#     use BesExternalArea BesExternalArea-* External
#     use GENBES * External (native_version=genbes-00-00-11)
#   use BesGenModule BesGenModule-* Generator
#   use BesServices BesServices-* Control
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
# use Babayaga Babayaga-* Generator
#   use BesPolicy BesPolicy-* 
#   use BesServices BesServices-* Control
#   use CERNLIB CERNLIB-* External
#   use GeneratorObject GeneratorObject-* Generator
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#   use BesCLHEP BesCLHEP-* External
# use RunEventNumberAlg RunEventNumberAlg-* Event
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-* External
#   use EventModel EventModel-* Event
# use Pi0EtaToGGRecAlg Pi0EtaToGGRecAlg-* Reconstruction
#   use BesPolicy BesPolicy-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use EventModel EventModel-* Event
#   use EvtRecEvent EvtRecEvent-* Event
#   use VertexFit VertexFit-* Analysis
# use DTagAlg DTagAlg-* Reconstruction
# use DummyLoadOldROOT DummyLoadOldROOT-* Utilities
#   use BesPolicy BesPolicy-01-* 
#   use BesROOT * External
#   use GaudiInterface GaudiInterface-01-* External
#
# Selection :
use CMT v1r25 (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib)
use efhlt efhlt-00-00-05 EventFilter/OnlineJointer (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Xt Xt-00-00-02 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesExternalArea BesExternalArea-00-00-22 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DimSetup DimSetup-00-00-06 DistBoss (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use XmlRpcC XmlRpcC-00-00-02 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GdmlManagement GdmlManagement-00-00-32 DetectorDescription (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use LCG_Platforms v1  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use LCG_Configuration v1  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use LCG_Settings v1  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use uuid v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use blas v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use lapack v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use cernlib v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use AIDA v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesAIDA BesAIDA-00-00-01 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use sqlite v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use mysql v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use MYSQL MYSQL-00-00-09 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use XercesC v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesFortranPolicy BesFortranPolicy-00-01-03  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CASTOR v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use CERNLIB CERNLIB-01-02-03 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HepPDT v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use HepMC v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use GENBES GENBES-00-00-08 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use v2llg v2llg-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use twogam twogam-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use tester tester-00-00-11 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use tauprd tauprd-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use sagerx sagerx-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use rhopi rhopi-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use radmu radmu-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use radgg radgg-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use radee radee-00-00-09 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ppgen ppgen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use p2mumu p2mumu-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use p2epem p2epem-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use p2bb p2bb-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use mugen mugen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use lund_crm lund_crm-00-00-09 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use lund lund-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use kstark kstark-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use kk2f kk2f-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use koralbe koralbe-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use howl howl-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use gamma2 gamma2-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use fsfgen fsfgen-00-00-09 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ffgen ffgen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use fff fff-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use epscat epscat-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use dssgen dssgen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use dsdgen dsdgen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ddprod ddprod-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ddgen ddgen-00-00-07 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use bhagen bhagen-00-00-08 External/GENBES (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CLHEP v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesCLHEP BesCLHEP-00-00-11 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesGeant4 BesGeant4-00-00-11 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GdmlToG4 GdmlToG4-00-00-11 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use xrootd v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use GCCXML v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use libunwind v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use tcmalloc v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use Python v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use pytools v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use Boost v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesBoost BesBoost-00-00-01 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ROOT v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesROOT BesROOT-00-00-08 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GdmlToRoot GdmlToRoot-00-00-14 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Reflex v1 LCG_Interfaces (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use GaudiPolicy v12r7  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use BesCxxPolicy BesCxxPolicy-00-01-01  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesPolicy BesPolicy-01-05-05  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GeneratorUtil GeneratorUtil-00-00-03 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ers ers-00-00-03 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use eformat eformat-00-00-04 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use IRawFile IRawFile-00-00-05 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DistBossUtil DistBossUtil-00-00-04 DistBoss (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ClientErrHandler ClientErrHandler-00-00-01 DistBoss (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RawFile RawFile-00-00-10 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use PhySim PhySim-00-00-10 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GenSim GenSim-00-00-07 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DecayChain DecayChain-00-00-03-slc6tag Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ProbTools ProbTools-00-00-01 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ProxyDict ProxyDict-00-00-01 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcRecoUtil MdcRecoUtil-00-01-08 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RootPolicy RootPolicy-00-01-03  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use AsciiDmp AsciiDmp-01-03-01 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use SimUtil SimUtil-00-00-37 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Identifier Identifier-00-02-15 DetectorDescription (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcGeneralClass EmcGeneralClass-00-00-04 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use facilities facilities-00-00-04 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use xmlBase xmlBase-00-00-03 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use rdbModel rdbModel-00-01-01 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GaudiKernel v28r8  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiUtils v4r6  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiAlg v14r6  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiPython v12r6  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiAud v9r9  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiCommonSvc v1r6  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiCoreSvc v1r5  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiSys v23r9  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use Gaudi v23r9  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use RootHistCnv v11r2  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use PartPropSvc v4r6  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiSvc v19r4  (/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiInterface GaudiInterface-01-03-07 External (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DummyLoadOldROOT DummyLoadOldROOT-00-00-01 Utilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use NetDataWriter NetDataWriter-00-00-03 DistBoss (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltStore HltStore-01-00-04 EventFilter/HltUtilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltDataTypes HltDataTypes-01-01-03 EventFilter/HltUtilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltTools HltTools-01-00-02 EventFilter/HltUtilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltSteerData HltSteerData-01-00-03 EventFilter/HltProcessor (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesKernel BesKernel-00-00-03 Control (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesServices BesServices-00-00-11 Control (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ZddEvent ZddEvent-00-00-04 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DataInfoSvc DataInfoSvc-00-00-03 Control (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcWaveform EmcWaveform-00-00-03 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofGeomSvc TofGeomSvc-00-00-12 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrigEvent TrigEvent-00-01-02 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesTimerSvc BesTimerSvc-00-00-12 Utilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use G4Geo G4Geo-00-00-12 DetectorDescription (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ROOTGeo ROOTGeo-00-00-15 DetectorDescription (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucGeomSvc MucGeomSvc-00-02-25 Muc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcRecGeoSvc EmcRecGeoSvc-01-01-07 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RelTable RelTable-00-00-02 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibData CalibData-00-01-09 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DatabaseSvc DatabaseSvc-00-00-24 Database (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use calibUtil calibUtil-00-00-39 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventModel EventModel-01-05-33 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RunEventNumberAlg RunEventNumberAlg-00-00-02 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltConfig HltConfig-01-01-05 EventFilter/HltProcessor (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrigMakerAlg TrigMakerAlg-00-00-02 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RawDataCnvBase RawDataCnvBase-01-00-03 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcTunningSvc MdcTunningSvc-00-00-26 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DedxCurSvc DedxCurSvc-00-00-17 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GeneratorObject GeneratorObject-00-01-05 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Babayaga Babayaga-00-00-25 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use GeneratorModule GeneratorModule-00-01-05 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Bhwide Bhwide-00-00-11 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesEvtGen BesEvtGen-00-01-96-slc6tag Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use KKMC KKMC-00-00-44-slc6tag Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TruSim TruSim-00-00-17 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use OfflineEventLoopMgr OfflineEventLoopMgr-00-00-15 Control (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofQCorrSvc TofQCorrSvc-00-00-10 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RawEvent RawEvent-00-03-19 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltEvent HltEvent-00-02-07 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucRawEvent MucRawEvent-00-02-02 Muc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcRawEvent EmcRawEvent-00-02-05 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofRawEvent TofRawEvent-00-02-07 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use LumiDigi LumiDigi-00-00-02 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcRawEvent MdcRawEvent-00-03-08 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MagneticField MagneticField-00-01-38  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use McTruth McTruth-00-02-19 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesGenModule BesGenModule-00-00-13 Generator (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use V2llg V2llg-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Twogam Twogam-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Tester Tester-00-01-04 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Tauprd Tauprd-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Sagerx Sagerx-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Rhopi Rhopi-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Radmu Radmu-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Radgg Radgg-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Radee Radee-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Ppgen Ppgen-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use P2mumu P2mumu-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use P2epem P2epem-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use P2bb P2bb-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Mugen Mugen-00-01-06 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Lund_crm Lund_crm-00-01-06 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Lund Lund-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Kstark Kstark-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use KK2f KK2f-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Koralbe Koralbe-00-01-06 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Howl Howl-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Gamma2 Gamma2-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Fsfgen Fsfgen-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Ffgen Ffgen-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Fff Fff-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Epscat Epscat-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Dssgen Dssgen-00-01-06 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Dsdgen Dsdgen-00-01-06 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Ddprod Ddprod-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Ddgen Ddgen-00-01-05 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Bhagen Bhagen-00-01-04 Generator/BesGenInterface (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcGeomSvc MdcGeomSvc-00-01-37 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcGeom MdcGeom-00-01-17 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BField BField-00-01-02 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DedxCorrecSvc DedxCorrecSvc-00-02-52 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EvTimeEvent EvTimeEvent-00-00-08 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use T0Dummy T0Dummy-00-00-05 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DstEvent DstEvent-00-02-50-p01 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibDataSvc CalibDataSvc-00-01-04 Calibration/CalibSvc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibXmlCnvSvc CalibXmlCnvSvc-00-01-01 Calibration/CalibSvc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibMySQLCnv CalibMySQLCnv-00-01-09 Calibration/CalibSvc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibTreeCnv CalibTreeCnv-00-01-12 Calibration/CalibSvc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibROOTCnv CalibROOTCnv-00-01-06 Calibration/CalibSvc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcCalibConstSvc EmcCalibConstSvc-00-00-10 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RealizationSvc RealizationSvc-00-00-38 Simulation/Realization (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use G4Svc G4Svc-00-01-51 Simulation (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcSim EmcSim-00-00-43 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibSvc CalibSvc-00-01-96 Calibration (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucCalibConstSvc MucCalibConstSvc-00-01-10 Muc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofSimSvc TofSimSvc-00-00-04 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EstTofCaliSvc EstTofCaliSvc-00-00-08 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MrpcCalibSvc MrpcCalibSvc-00-00-02 Mrpc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcCalibFunSvc MdcCalibFunSvc-00-03-16 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcSim MdcSim-00-00-70 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofQElecSvc TofQElecSvc-00-00-05 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofCaliSvc TofCaliSvc-00-01-05 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofSim TofSim-00-02-12 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ExtEvent ExtEvent-00-00-13 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ReconEvent ReconEvent-00-00-04 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucRecEvent MucRecEvent-00-02-51 Muc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RootEventData RootEventData-00-03-80 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesVisLib BesVisLib-00-04-04 EventDisplay (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesVisClient BesVisClient-00-04-04 EventDisplay (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcRecEventModel EmcRecEventModel-01-01-18 Emc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcTimeRec EmcTimeRec-00-00-03 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofRecEvent TofRecEvent-00-02-14 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcRecEvent MdcRecEvent-00-05-14 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DstMakerAlg DstMakerAlg-00-02-28 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltMakerAlg HltMakerAlg-00-00-06 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventNavigator EventNavigator-00-01-03 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucCalibAlg MucCalibAlg-00-02-16 Muc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucSim MucSim-00-01-03 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesSim BesSim-00-01-22 Simulation/BOOST (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MucRecAlg MucRecAlg-00-03-08 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrkBase TrkBase-00-01-12 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrkFitter TrkFitter-00-01-11 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcData MdcData-00-01-27 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcUtilitySvc MdcUtilitySvc-00-00-07 Mdc/MdcCheckUtil (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrackUtil TrackUtil-00-00-08 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrkExtAlg TrkExtAlg-00-00-60.bak Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcTables MdcTables-00-00-11 Mdc (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RawDataProviderSvc RawDataProviderSvc-00-03-40 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltAlgorithms HltAlgorithms-01-03-01 EventFilter/HltUtilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use HltSteering HltSteering-01-01-05 EventFilter/HltProcessor (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RawDataCnv RawDataCnv-00-04-29-tag Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EFServices EFServices-00-01-09 EventFilter/OnlineJointer (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ESController ESController-00-01-05 EventFilter/OnlineJointer (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventFilter EventFilter-02-00-10  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TrkReco TrkReco-00-08-59-patch4-slc6tag Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcFastTrkAlg MdcFastTrkAlg-00-04-09 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EmcRec EmcRec-01-02-52 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofRec TofRec-00-04-30 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use TofEnergyRec TofEnergyRec-00-00-08 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use tofcalgsec tofcalgsec-00-02-08 Tof (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EvtRecEvent EvtRecEvent-00-01-05 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Trigger Trigger-00-01-05 Trigger (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RootCnvSvc RootCnvSvc-02-01-05 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesVisAlg BesVisAlg-00-01-07 EventDisplay (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventWriter EventWriter-00-00-18 EvtPreSelect (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesEventMixer BesEventMixer-00-00-34 Simulation (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RootIO RootIO-00-01-31 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use SimplePIDSvc SimplePIDSvc-00-00-11 Utilities (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventAssembly EventAssembly-00-00-17 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcPrintSvc MdcPrintSvc-00-00-01 Mdc/MdcCheckUtil (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcTrkRecon MdcTrkRecon-00-03-45 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcxReco MdcxReco-00-01-59 Reconstruction/MdcPatRec (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EsTimeAlg EsTimeAlg-00-02-52 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use MdcDedxAlg MdcDedxAlg-00-06-57 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Reconstruction Reconstruction-00-00-20  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use ParticleID ParticleID-00-04-55 Analysis (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use VertexFit VertexFit-00-02-74 Analysis (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Pi0EtaToGGRecAlg Pi0EtaToGGRecAlg-00-00-10 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use Gam4pikpAlg Gam4pikpAlg-00-00-08 Analysis/ControlSamples (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DimuPreSelect DimuPreSelect-00-00-09 EvtPreSelect (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use EventPreSelect EventPreSelect-00-00-21 EvtPreSelect (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DiGamAlg DiGamAlg-00-00-24 DQA (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use LumTauAlg LumTauAlg-00-00-04  (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BesDChain BesDChain-00-00-14 Event (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DTagAlg DTagAlg-00-00-51 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use CalibEventSelect CalibEventSelect-00-00-15 EvtPreSelect (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use DTagSkim DTagSkim-00-00-06 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use BeamParamsAlg BeamParamsAlg-00-00-10 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use VeeVertexAlg VeeVertexAlg-00-00-05 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use PrimaryVertexAlg PrimaryVertexAlg-00-00-04 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use KalFitAlg KalFitAlg-00-07-55-p03 Reconstruction (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use AbsCor AbsCor-00-00-28 Analysis/PhotonCor (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use PipiJpsiAlg PipiJpsiAlg-00-00-03 Analysis/Physics/PsiPrime (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
use RhopiAlg RhopiAlg-00-00-25 Analysis/Physics (/afs/ihep.ac.cn/users/x/xiuql/scratchfs/bes/jpsi2invi/v0.1)
use ExHelloWorld ExHelloWorld-00-00-03 BesExamples (/afs/ihep.ac.cn/bes3/offline/Boss/6.6.5)
----------> tags
CMTv1 (from CMTVERSION)
CMTr25 (from CMTVERSION)
CMTp0 (from CMTVERSION)
Linux (from uname) package [CMT LCG_Platforms BesPolicy] implies [Unix host-linux]
x86_64-linux (from HOSTTYPE)
x86_64-slc6-gcc46-opt (from CMTCONFIG) package [LCG_Platforms] implies [target-x86_64 target-slc6 target-gcc46 target-opt]
jpsi2invi_no_config (from PROJECT) excludes [jpsi2invi_config]
jpsi2invi_root (from PROJECT) excludes [jpsi2invi_no_root]
jpsi2invi_cleanup (from PROJECT) excludes [jpsi2invi_no_cleanup]
jpsi2invi_scripts (from PROJECT) excludes [jpsi2invi_no_scripts]
jpsi2invi_prototypes (from PROJECT) excludes [jpsi2invi_no_prototypes]
jpsi2invi_with_installarea (from PROJECT) excludes [jpsi2invi_without_installarea]
jpsi2invi_with_version_directory (from PROJECT) excludes [jpsi2invi_without_version_directory]
jpsi2invi (from PROJECT)
BOSS_no_config (from PROJECT) excludes [BOSS_config]
BOSS_root (from PROJECT) excludes [BOSS_no_root]
BOSS_cleanup (from PROJECT) excludes [BOSS_no_cleanup]
BOSS_scripts (from PROJECT) excludes [BOSS_no_scripts]
BOSS_no_prototypes (from PROJECT) excludes [BOSS_prototypes]
BOSS_with_installarea (from PROJECT) excludes [BOSS_without_installarea]
BOSS_with_version_directory (from PROJECT) excludes [BOSS_without_version_directory]
GAUDI_no_config (from PROJECT) excludes [GAUDI_config]
GAUDI_root (from PROJECT) excludes [GAUDI_no_root]
GAUDI_cleanup (from PROJECT) excludes [GAUDI_no_cleanup]
GAUDI_scripts (from PROJECT) excludes [GAUDI_no_scripts]
GAUDI_prototypes (from PROJECT) excludes [GAUDI_no_prototypes]
GAUDI_with_installarea (from PROJECT) excludes [GAUDI_without_installarea]
GAUDI_without_version_directory (from PROJECT) excludes [GAUDI_with_version_directory]
LCGCMT_no_config (from PROJECT) excludes [LCGCMT_config]
LCGCMT_no_root (from PROJECT) excludes [LCGCMT_root]
LCGCMT_cleanup (from PROJECT) excludes [LCGCMT_no_cleanup]
LCGCMT_scripts (from PROJECT) excludes [LCGCMT_no_scripts]
LCGCMT_prototypes (from PROJECT) excludes [LCGCMT_no_prototypes]
LCGCMT_without_installarea (from PROJECT) excludes [LCGCMT_with_installarea]
LCGCMT_with_version_directory (from PROJECT) excludes [LCGCMT_without_version_directory]
x86_64 (from package CMT) package [LCG_Platforms] implies [host-x86_64]
sl65 (from package CMT)
gcc463 (from package CMT)
Unix (from package CMT) package [LCG_Platforms] implies [host-unix] excludes [WIN32 Win32]
c_native_dependencies (from package CMT) activated GaudiPolicy
cpp_native_dependencies (from package CMT) activated GaudiPolicy
target-unix (from package LCG_Settings) activated LCG_Platforms
target-x86_64 (from package LCG_Settings) activated LCG_Platforms
target-gcc46 (from package LCG_Settings) package [LCG_Platforms] implies [target-gcc4 target-lcg-compiler lcg-compiler] activated LCG_Platforms
host-x86_64 (from package LCG_Settings) activated LCG_Platforms
target-slc (from package LCG_Settings) package [LCG_Platforms] implies [target-linux] activated LCG_Platforms
target-gcc (from package LCG_Settings) activated LCG_Platforms
target-gcc4 (from package LCG_Settings) package [LCG_Platforms] implies [target-gcc] activated LCG_Platforms
target-lcg-compiler (from package LCG_Settings) activated LCG_Platforms
host-linux (from package LCG_Platforms) package [LCG_Platforms] implies [host-unix]
host-unix (from package LCG_Platforms)
target-opt (from package LCG_Platforms)
target-slc6 (from package LCG_Platforms) package [LCG_Platforms] implies [target-slc]
target-linux (from package LCG_Platforms) package [LCG_Platforms] implies [target-unix]
lcg-compiler (from package LCG_Platforms)
ROOT_GE_5_15 (from package LCG_Configuration)
ROOT_GE_5_19 (from package LCG_Configuration)
HasAthenaRunTime (from package BesPolicy)
ROOTBasicLibs (from package BesROOT)
----------> CMTPATH
# Add path /afs/ihep.ac.cn/users/x/xiuql/scratchfs/bes/jpsi2invi/v0.1 from initialization
# Add path /afs/ihep.ac.cn/bes3/offline/Boss/6.6.5 from initialization
# Add path /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9 from initialization
# Add path /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a from initialization
