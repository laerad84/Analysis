{
  
  gSystem -> Load("libE14Data.so");

  E14Fill a;
  a.SetDebugMode();
  a.SetInputDirectory("../conv");
  a.SetRunNumber(1512);
  a.SetOutputFile("out.root");
  a.FillFromConv();

}
