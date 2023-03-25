using System.Globalization;
NumberFormatInfo setPrecision = new NumberFormatInfo();    
setPrecision.NumberDecimalDigits = 8;
// See https://aka.ms/new-console-template for more information
double[] X1(double t) { return new [] { Math.Cos(t), Math.Sin(t) }; }
double[] X2(double t) { return new [] { 2 * Math.Cos(t), 2 * Math.Sin(t) }; }
double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), Math.Cos(t) }; }
double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), 2 * Math.Cos(t) }; }
double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -1 * Math.Sin(t) }; }
double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -2 * Math.Sin(t) }; }
// double[] X1(double t) { return new [] { Math.Cos(t), 0.5 * Math.Sin(t) }; }
// double[] X2(double t) { return new [] { 2 * Math.Cos(t), Math.Sin(t) }; }
// double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), 0.5 * Math.Cos(t) }; }
// double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), Math.Cos(t) }; }
// double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -0.5 * Math.Sin(t) }; }
// double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -1 * Math.Sin(t) }; }
double Eps = 0.000001;

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

double K11(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (Der2X1(t1)[0] * Der1X1(t1)[1] - Der2X1(t1)[1] * Der1X1(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2));
    }

    double kernelNumerator = (X1(t2)[0] - X1(t1)[0]) * VGamma1(t1)[0] + (X1(t2)[1] - X1(t1)[1]) * VGamma1(t1)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0) /
           kernelDenominator;
}

double K12(double t1, double t2)
{
    double kernelLeftValue = (VGamma1(t1)[0] * VGamma2(t2)[0] + VGamma1(t1)[1] * VGamma2(t2)[1]) /
                             Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 2);
    double kernelRightValueNominator = ((X1(t1)[0] - X2(t2)[0]) * VGamma2(t2)[0] + (X1(t1)[1] - X2(t2)[1]) * VGamma2(t2)[1]) *
                                       ((X1(t1)[0] - X2(t2)[0]) * VGamma1(t1)[0] + (X1(t1)[1] - X2(t2)[1]) * VGamma1(t1)[1]);
    double kernelRightValue = kernelRightValueNominator / 
                              Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 4);
    
    return (kernelLeftValue - 2 * kernelRightValue) * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0);
}

double K21(double t1, double t2)
{
    double kernelLeftValue = Math.Log(1 / GetEuclideanDistance(X2(t1)[0], X1(t2)[0], X2(t1)[1], X1(t2)[1]));
    double kernelRightValue = GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0);
    return kernelLeftValue * kernelRightValue;
}

double K22(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (-1 * Der1X2(t1)[0] * Der2X2(t1)[1] + Der1X2(t1)[1] * Der2X2(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X2(t1)[0], 0, Der1X2(t1)[1], 0), 2));
    }

    double kernelNumerator = (X2(t1)[0] - X2(t2)[0]) * VGamma2(t2)[0] + (X2(t1)[1] - X2(t2)[1]) * VGamma2(t2)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X2(t1)[0], X2(t2)[0], X2(t1)[1], X2(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
        kernelDenominator;
}

double H11(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (-1 * Der1X1(t1)[0] * Der2X1(t1)[1] + Der1X1(t1)[1] * Der2X1(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2));
    }

    double kernelNumerator = (X1(t1)[0] - X1(t2)[0]) * VGamma1(t2)[0] + (X1(t1)[1] - X1(t2)[1]) * VGamma1(t2)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0) /
           kernelDenominator;
}

double H12(double t1, double t2)
{
    double kernelLeftValue = Math.Log(1 / GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]));
    double kernelRightValue = GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0);
    return kernelLeftValue * kernelRightValue;
}

double H21(double t1, double t2)
{
    double kernelLeftValue = (VGamma2(t1)[0] * VGamma1(t2)[0] + VGamma2(t1)[1] * VGamma1(t2)[1]) /
                             Math.Pow(GetEuclideanDistance(X2(t1)[0], X1(t2)[0], X2(t1)[1], X1(t2)[1]), 2);
    double kernelRightValueNominator = ((X2(t1)[0] - X1(t2)[0]) * VGamma1(t2)[0] + (X2(t1)[1] - X1(t2)[1]) * VGamma1(t2)[1]) *
                                       ((X2(t1)[0] - X1(t2)[0]) * VGamma2(t1)[0] + (X2(t1)[1] - X1(t2)[1]) * VGamma2(t1)[1]);
    double kernelRightValue = kernelRightValueNominator / 
                              Math.Pow(GetEuclideanDistance(X2(t1)[0], X1(t2)[0], X2(t1)[1], X1(t2)[1]), 4);
    
    return (kernelLeftValue - 2 * kernelRightValue) * GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0); 
}

double H22(double t1, double t2)
{
    if (Math.Abs(t1 - t2) < Eps)
    {
        return (Der2X2(t1)[0] * Der1X2(t1)[1] - Der2X2(t1)[1] * Der1X2(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X2(t1)[0], 0, Der1X2(t1)[1], 0), 2));
    }

    double kernelNumerator = (X2(t2)[0] - X2(t1)[0]) * VGamma2(t1)[0] + (X2(t2)[1] - X2(t1)[1]) * VGamma2(t1)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X2(t1)[0], X2(t2)[0], X2(t1)[1], X2(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
           kernelDenominator;
}


double GetEuclideanDistance(double x21, double x11, double x22, double x12)
{
    return Math.Sqrt(Math.Pow(x21 - x11, 2) + Math.Pow(x22 - x12, 2));
}

double Rj(double t, double tj, double n)
{
    double result = 0;
    
    for (int m = 1; m < n; m++)
    {
        result += 1.0 / m * Math.Cos(m * (t - tj)) + 1 / (2 * n) * Math.Cos(n * (t - tj));

    }

    return - result / n;
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
        // s1 = s1 + uValues[i] * Math.Log(1 / GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1])) * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0);
        // s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0) *
        //       (((x1 - X2(t)[0]) * VGamma2(t)[0] + (x2 - X2(t)[1]) * VGamma2(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1]), 2));
        s1 = s1 + uValues[i] * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0) *
            (((x1 - X1(t)[0]) * VGamma1(t)[0] + (x2 - X1(t)[1]) * VGamma1(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1]), 2));
        s2 = s2 + uValues[i+2*N] * Math.Log(1 / GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1])) * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0);

    }

    return s1 / (2 * N) + s2 / (2 * N);
}

double[] GetApproximatedUOnGamma1(double[] uValues, int N)
{
    double[] result = new double[2*N];
    
    for (int k = 0; k < 2*N; k++)
    {
        double tk = k * Math.PI / N;
        double s1 = 0, s2 = 0;
        for (int i = 0; i < 2*N; i++)
        {
            double ti = i * Math.PI / N;
            s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(ti)[0], 0, Der1X2(ti)[1], 0) *
                ((X1(tk)[0] - X2(ti)[0]) * VGamma2(ti)[0] + (X1(tk)[1] - X2(ti)[1]) * VGamma2(ti)[1]) / Math.Pow(GetEuclideanDistance(X1(tk)[0], X2(ti)[0], X1(tk)[1], X2(ti)[1]), 2);

            double l = 0;
            if (Math.Abs(ti - tk) < Eps)
            {
                l = Math.Log(Math.Pow(GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0), 2));
            }
            else
            {
                l = Math.Log(Math.Pow(GetEuclideanDistance(X1(tk)[0], X1(ti)[0], X1(tk)[1], X1(ti)[1]), 2) / (4 * Math.Pow(Math.Sin((tk - ti)/2), 2)));
            }

            s1 = s1 + uValues[i] * GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0) * (-l/(4*N) - Rj(ti,tk, N)/2);
        }

        result[k] = s1 + s2 / (2 * N);
    }


    return result;
}

void Solve(int N1)
{
    int N = N1;
    double[] H_F_Values = new double [4*N];
    double[,] kernelMatrix = new double [4*N, 4*N];
    // for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? 0 : 1; }
    for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? 1 : 0; }

    for (int i = 0; i < 2*N; i++)
    {
        double ti = i * Math.PI / N;
        for (int j = 0; j < 2*N; j++)
        {
            double tj = j * Math.PI / N;
            // kernelMatrix[i, j] = K11(ti, tj)/(2*N);
            // kernelMatrix[i, 2*N + j] = K12(ti, tj)/(2*N);
            // kernelMatrix[2*N + i, j] = K21(ti, tj)/(2*N);
            // kernelMatrix[2*N + i, 2*N + j] = K22(ti, tj)/(2*N);
            kernelMatrix[i, j] = H11(ti, tj)/(2*N);
            kernelMatrix[i, 2*N + j] = H12(ti, tj)/(2*N);
            kernelMatrix[2*N + i, j] = H21(ti, tj)/(2*N);
            kernelMatrix[2*N + i, 2*N + j] = H22(ti, tj)/(2*N);
        }

        // kernelMatrix[i, i] = kernelMatrix[i, i] - 0.5;
        // kernelMatrix[2*N + i, 2*N + i] = kernelMatrix[2*N + i, 2*N + i] - 0.5;
        kernelMatrix[i, i] = kernelMatrix[i, i] + 0.5;
        kernelMatrix[2*N + i, 2*N + i] = kernelMatrix[2*N + i, 2*N + i] + 0.5;
    }

    double[,] kernelMatrixExtended = new double[4*N, 4*N + 1];
    for (int i = 0; i < 4*N; i++) { kernelMatrixExtended[i, 4*N] = H_F_Values[i]; }
    for (int i = 0; i < 4*N; i++) { for (int j = 0; j < 4*N; j++) { kernelMatrixExtended[i, j] = kernelMatrix[i, j]; } }
    
    // for (int i = 0; i < 4 * N; i++) { for (int j = 0; j < 4 * N + 1; j++) { Console.Write( kernelMatrixExtended[i,j].ToString("N", setPrecision) + " "); } Console.WriteLine(); } Console.WriteLine();
    // double[,] testCase = new double[3, 4] { { 1, 9, -5, -32 }, { -3, -5, -5, -10 }, { -2, -7, 1, 13 } }; double[] ans1 = Gauss(testCase, 3); for (int j = 0; j < 3; j++) { Console.Write(ans1[j] + " "); }
    
    double[] ans = Gauss(kernelMatrixExtended, 4 * N);
    // Console.WriteLine("Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(ans[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    Console.WriteLine("~U(1.5, 0) Value: " + GetApproximatedU(1.5, 0, ans, N).ToString("N", setPrecision) + " ");
    Console.WriteLine("~U(0, 0.75) Value: " + GetApproximatedU(0, 0.75, ans, N).ToString("N", setPrecision) + " ");


    // double[] uOnGamma1 = GetApproximatedUOnGamma1(ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(uOnGamma1[j].ToString("N", setPrecision) + " "); }
}

Solve(4);
Console.WriteLine();
Solve(8);
Console.WriteLine();
Solve(16);
Console.WriteLine();
Solve(32);

Console.WriteLine("\nEND");