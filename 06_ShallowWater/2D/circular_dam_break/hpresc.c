double hpresc(double x, double y)
{
    return ((x-20.0)*(x-20.0) + (y-20.0)*(y-20.0)) <= 6.25 ? 12.5 : 0.5;
}