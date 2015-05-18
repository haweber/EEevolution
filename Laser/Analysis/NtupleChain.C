{
  TChain data("x");

  //all raw DATA files are added to chain so that one runs over all
  //TString path = "/shome/donega/ECAL/DATA";
  Int_t n = 0;
  n+=data.Add("/shome/haweber/ECAL/DataLaser/Antu-data_00173853-00173921.root",0);
  n+=data.Add("/shome/haweber/ECAL/DataLaser/Antu-data_00174038-00175994.root",0);
  n+=data.Add("/shome/haweber/ECAL/DataLaser/Antu-data_00176002-00176982.root",0);
  n+=data.Add("/shome/haweber/ECAL/DataLaser/Antu-data_00177026-00177954.root",0);

  cout << "#files = " << n << endl;


}