// See https://aka.ms/new-console-template for more information
double[] X1(double t) { return new [] { Math.Cos(t), Math.Sin(t) }; }
double[] X2(double t) { return new [] { 2 * Math.Cos(t), 2 * Math.Sin(t) }; }
double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), Math.Cos(t) }; }
double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), 2 * Math.Cos(t) }; }
double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -1 * Math.Sin(t) }; }
double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -2 * Math.Sin(t) }; }
double Eps = 0.00000001;

double[] VGamma2(double t)
{
    return new []
    {
        Der1X2(t)[1] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0),
        -1 * Der1X2(t)[0] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0)
    };
}
double[] VGamma1(double t)
{
    return new []
    {
        Der1X1(t)[1] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0),
        -1 * Der1X1(t)[0] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0)
    };
}

double Kernel11(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (-1 * Der1X1(t1)[0] * Der2X1(t1)[1] + Der1X1(t1)[1] * Der2X1(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 3));
    }

    double kernelNumerator = (-X1(t1)[0] + X1(t2)[0]) * VGamma1(t1)[0] + (-X1(t1)[1] + X1(t2)[1]) * VGamma1(t1)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0) /
           kernelDenominator;
}

double Kernel12(double t1, double t2)
{
    double kernelLeftValue = (VGamma2(t2)[0] * VGamma1(t1)[0] + VGamma2(t2)[1] * VGamma1(t1)[1]) /
                             Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 2);
    double kernelRightValueNominator = ((X1(t1)[0] - X2(t2)[0]) * VGamma2(t2)[0] + (X1(t1)[1] - X2(t2)[1]) * VGamma2(t2)[1]) *
                                       ((X1(t1)[0] - X2(t2)[0]) * VGamma1(t1)[0] + (X1(t1)[1] - X2(t2)[1]) * VGamma1(t1)[1]);
    double kernelRightValue = kernelRightValueNominator / 
                              Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 4);
    
    return (kernelLeftValue - 2 * kernelRightValue) * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0);
}

double Kernel21(double t1, double t2)
{
    double kernelLeftValue = Math.Log(1 / GetEuclideanDistance(X2(t1)[0], X1(t2)[0], X2(t1)[1], X1(t2)[1]));
    double kernelRightValue = GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0);
    return kernelLeftValue * kernelRightValue;
}

double Kernel22(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (-1 * Der1X2(t1)[0] * Der2X2(t1)[1] + Der1X2(t1)[1] * Der2X2(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X2(t1)[0], 0, Der1X2(t1)[1], 0), 3));
    }

    double kernelNumerator = (X2(t1)[0] - X2(t2)[0]) * VGamma2(t2)[0] + (X2(t1)[1] - X2(t2)[1]) * VGamma2(t2)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X2(t1)[0], X2(t2)[0], X2(t1)[1], X2(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
        kernelDenominator;
}

double GetEuclideanDistance(double x21, double x11, double x22, double x12)
{
    return Math.Sqrt(Math.Pow(x21 - x11, 2) + Math.Pow(x22 - x12, 2));
}

double[] Gauss (double[,] a, int N) {
    int n = N;
    int m = N;
 
    int[] where = new int[m];
    for (int i=0; i<m; i++) { where[i] = -1; }
    
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (Math.Abs(a[i, col]) > Math.Abs(a[sel, col]))
                sel = i;
        if (Math.Abs (a[sel, col]) < Eps)
            continue;
        for (int i = col; i <= m; ++i)
            (a[sel, i], a[row, i]) = (a[row, i], a[sel, i]);
        where[col] = row;
 
        for (int i=0; i<n; ++i)
            if (i != row) {
                double c = a[i, col] / a[row, col];
                for (int j=col; j<=m; ++j)
                    a[i, j] -= a[row, j] * c;
            }
        ++row;
    }
    
    // for (int i = 0; i < m; i++) { Console.Write(where[i] + " "); } Console.WriteLine();
    
    double[] ans = new double[m];
    for (int i=0; i<m; i++) { ans[i] = 0; }
    
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i], m] / a[where[i], i];
    for (int i=0; i<n; ++i) {
        double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i, j];
        if (Math.Abs(sum - a[i, m]) > Eps)
            return new double[] { };;
    }
 
    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return new double[] { };
    return ans;
}

double GetApproximatedU(double x1, double x2, double[] uValues, int N)
{
    double s1 = 0, s2 = 0;
    for (int i = 0; i < 2*N; i++)
    {
        double t = i * Math.PI / N;
        s1 = s1 + uValues[i] * Math.Log(1 / GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1])) * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0);
        s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0) *
              (((x1 - X2(t)[0]) * VGamma2(t)[0] + (x2 - X2(t)[1]) * VGamma2(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1]), 2));
    }

    return s1 / (2 * N) + s2 / (2 * N);
}

void Solve(int N1)
{
    int N = N1;
    double[] H_F_Values = new double [4*N];
    double[,] kernelMatrix = new double [4*N, 4*N];
    for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? 0 : 1; }

    for (int i = 0; i < 2*N; i++)
    {
        for (int j = 0; j < 2*N; j++)
        {
            kernelMatrix[i, j] = Kernel11(i * Math.PI/(2*N), j * Math.PI/(2*N))/(2*N);
            kernelMatrix[i, 2*N + j] = Kernel12(i * Math.PI/(2*N), 2*N*j * Math.PI/(2*N))/(2*N);
            kernelMatrix[2*N + i, j] = Kernel21(2*N*i * Math.PI/(2*N), j * Math.PI/(2*N))/(2*N);
            kernelMatrix[2*N + i, 2*N + j] = Kernel22(2*N*i * Math.PI/(2*N), 2*N*j * Math.PI/(2*N))/(2*N);
        }

        kernelMatrix[i, i] = kernelMatrix[i, i] - 0.5;
        kernelMatrix[2*N + i, 2*N + i] = kernelMatrix[2*N + i, 2*N + i] - 0.5;
    }

    double[,] kernelMatrixExtended = new double[4*N, 4*N + 1];
    for (int i = 0; i < 4*N; i++) { kernelMatrixExtended[i, 4*N] = H_F_Values[i]; }
    for (int i = 0; i < 4*N; i++) { for (int j = 0; j < 4*N; j++) { kernelMatrixExtended[i, j] = kernelMatrix[i, j]; } }
    
    // for (int i = 0; i < 4 * N; i++) { for (int j = 0; j < 4 * N + 1; j++) { Console.Write(kernelMatrixExtended[i,j] + " "); } Console.WriteLine(); } Console.WriteLine();
    // double[,] testCase = new double[2, 3] { { 1, 2, 3 }, { 4, 5, 6 } }; double[] ans = Gauss(testCase, 2); for (int j = 0; j < 2; j++) { Console.Write(ans[j] + " "); }
    
    double[] ans = Gauss(kernelMatrixExtended, 4 * N);
    Console.WriteLine("Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(ans[j] + " "); }
    Console.WriteLine("\n");
    
    Console.WriteLine("~U Value: " + GetApproximatedU(1.5, 1.5, ans, N));
}

Solve(4);
Console.WriteLine();
Solve(8);
Console.WriteLine();
Solve(16);
Console.WriteLine();
Solve(32);

Console.WriteLine("\nEND");