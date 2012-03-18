// gcc -shared -fPIC 1.c -o lib1.so
double
fun(double(*callback)(double))
{
  return callback(3);
}
