void TestMatrix(){
  TMatrixD T(5,5); 
  TArrayD data(25);
  for( int i = 0; i< 25; i++){
    const int ir = i/5;
    const int ic = i%5;
    data[i] = 1./(ir+ic);
  }
  T.SetMatrixArray( data.GetArray() );
  T.Invert();
  T.Print();
}
