double zb(double x, double y)
{
    return (8.0 < x || x < 12.0) ? 0.2 - 0.05*(x-10)*(x-10) : 0.0;
}