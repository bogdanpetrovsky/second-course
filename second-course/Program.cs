using System.Globalization;
using second_course;

NumberFormatInfo setPrecision = new NumberFormatInfo();    
setPrecision.NumberDecimalDigits = 8;
// See https://aka.ms/new-console-template for more information

void Solve()
{
    int N = 32;
    Console.WriteLine("N:" + N);
    double[] exactSolution = new double [N*2];
    
    double[] f1Values = new double [N*2];
    double[] f2Values = new double [N*2];
    double[] g1Values = new double [N*2];
    double[] g2Values = new double [N*2];

    for (int i = 0; i < N * 2; i++)
    {
        double ti = i * Math.PI / N;
        exactSolution[i] = FunctionHelper.F1(ti);
        f1Values[i] = FunctionHelper.F1(ti);
        f2Values[i] = FunctionHelper.F1(ti);
        g1Values[i] = FunctionHelper.G2(ti);
        g2Values[i] = FunctionHelper.G2(ti);
    }
    
    // double[,] testCase = new double[3, 4] { { 1, 9, -5, -32 }, { -3, -5, -5, -10 }, { -2, -7, 1, 13 } }; double[] ans1 = FunctionHelper.Gauss(testCase, 3); for (int j = 0; j < 3; j++) { Console.Write(ans1[j] + " "); }
    ND_Solver ndSolver = new ND_Solver();
    DN_Solver dnSolver = new DN_Solver();
    
    double solutionError = Double.MaxValue;

    while (solutionError > FunctionHelper.ErrorEps)
    {
        double[] DN_Ans = dnSolver.Solve(N, f1Values, g2Values);
        double[] DN_derUOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(DN_Ans, N);
        g1Values = DN_derUOnGamma1;
        
        double[] ND_Ans = ndSolver.Solve(N, f2Values, g1Values);
        double[] ND_UOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(ND_Ans, N);
        f1Values = ND_UOnGamma1;

        solutionError = FunctionHelper.GetMaxError(f1Values, exactSolution);
        Console.WriteLine(solutionError);
    }
    
    
    
    // Console.WriteLine("ND Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(ND_Ans[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    // Console.WriteLine("DN Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(DN_Ans[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    // Console.WriteLine("ND ~U(1.5, 0) Value: " + FunctionHelper.GetApproximatedU(1.5, 0, ND_Ans, N, false).ToString("N", setPrecision) + " ");
    // Console.WriteLine("DN ~U(1.5, 0) Value: " + FunctionHelper.GetApproximatedU(1.5, 0, DN_Ans, N, true).ToString("N", setPrecision) + " ");

    // Console.WriteLine("~U(0, 1.25) Value: " + FunctionHelper.GetApproximatedU(0, 1.25, DN_Ans, N, true).ToString("N", setPrecision) + " ");

    // double[] ND_uOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(ND_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(ND_uOnGamma1[j].ToString("N", setPrecision) + " "); }
    
    // double[] DN_uOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(DN_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(DN_uOnGamma1[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    // double[] DN_deruOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(DN_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(Math.Abs(DN_deruOnGamma1[j]-FunctionHelper.du_G1(j*Math.PI / N)).ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    // double[] ND_deruOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(ND_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(Math.Abs(ND_deruOnGamma1[j]-FunctionHelper.du_G1(j*Math.PI / N)).ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
}

Solve();
Console.WriteLine();
// Solve(8);
// Console.WriteLine();
// Solve(16);
// Console.WriteLine();
// Solve(32);
// Console.WriteLine();
// Solve(64);
// Console.WriteLine();

Console.WriteLine("\nEND");

public static class FunctionHelper
{
    // double[] X1(double t) { return new [] { Math.Cos(t), Math.Sin(t) }; }
    // double[] X2(double t) { return new [] { 2 * Math.Cos(t), 2 * Math.Sin(t) }; }
    // double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), Math.Cos(t) }; }
    // double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), 2 * Math.Cos(t) }; }
    // double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -1 * Math.Sin(t) }; }
    // double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -2 * Math.Sin(t) }; }
    // double[] Der3X1(double t) { return new [] { Math.Sin(t), -1 * Math.Cos(t) }; }
    // double[] Der3X2(double t) { return new [] { 2 * Math.Sin(t), -2 * Math.Cos(t) }; }
    public static double[] X1(double t) { return new [] { Math.Cos(t), 0.5 * Math.Sin(t) }; }
    public static double[] X2(double t) { return new [] { 1.5 * Math.Cos(t), Math.Sin(t) }; }
    public static double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), 0.5 * Math.Cos(t) }; }
    public static double[] Der1X2(double t) { return new [] { -1.5 * Math.Sin(t), Math.Cos(t) }; }
    public static double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -0.5 * Math.Sin(t) }; }
    public static double[] Der2X2(double t) { return new [] { -1.5 * Math.Cos(t), -1 * Math.Sin(t) }; }
    public static double[] Der3X1(double t) { return new [] { Math.Sin(t), -0.5 * Math.Cos(t) }; }
    public static double[] Der3X2(double t) { return new [] { 1.5 * Math.Sin(t), -1 * Math.Cos(t) }; }
    public static double Eps = 0.000001;
    public static double ErrorEps = 0.01;
    
    public static double[] VGamma2(double t)
    {
        return new []
        {
            Der1X2(t)[1] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0),
            -1 * Der1X2(t)[0] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0)
        };
    }
    public static double[] VGamma1(double t)
    {
        return new []
        {
            Der1X1(t)[1] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0),
            -1 * Der1X1(t)[0] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0)
        };
    }
    
    public static double GetEuclideanDistance(double x21, double x11, double x22, double x12)
    {
        return Math.Sqrt(Math.Pow(x21 - x11, 2) + Math.Pow(x22 - x12, 2));
    }
    
    public static double[] Gauss (double[,] a, int N) {
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
                return ans;
        }
 
        for (int i=0; i<m; ++i)
            if (where[i] == -1)
                return ans;
        return ans;
    }
    
    public static double Rj(double t, double tj, double n)
    {
        double result = 0;
        
        for (int m = 1; m < n; m++)
        {
            result += 1.0 / m * Math.Cos(m * (t - tj)) + 1 / (2 * n) * Math.Cos(n * (t - tj));

        }

        return - result / n;
    }

    public static double GetApproximatedU(double x1, double x2, double[] uValues, int N, bool isReversedProblem)
    {
        double s1 = 0, s2 = 0;
        for (int i = 0; i < 2*N; i++)
        {
            double t = i * Math.PI / N;
            if (isReversedProblem)
            {
                s1 = s1 + uValues[i] * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0) *
                    (((x1 - X1(t)[0]) * VGamma1(t)[0] + (x2 - X1(t)[1]) * VGamma1(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1]), 2));
                s2 = s2 + uValues[i+2*N] * Math.Log(1 / GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1])) * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0);
            }
            else
            {
                s1 = s1 + uValues[i] * Math.Log(1 / GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1])) * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0);
                s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0) *
                      (((x1 - X2(t)[0]) * VGamma2(t)[0] + (x2 - X2(t)[1]) * VGamma2(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1]), 2));
            }
        }

        return s1 / (2 * N) + s2 / (2 * N);
    }

    public static double[] GetApproximatedUOnGamma1(double[] uValues, int N)
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

    public static double HWaved(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < Eps)
        {
            return -1.0 / 6
                   +
                   (Der1X1(t1)[0] * Der3X1(t1)[0] + Der1X1(t2)[1] * Der3X1(t1)[1]) /
                   (3 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2))
                   +
                   Math.Pow(GetEuclideanDistance(Der2X1(t1)[0], 0, Der2X1(t1)[1], 0), 2) /
                   (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2))
                   -
                   Math.Pow(Der1X1(t1)[0] * Der2X1(t1)[0] + Der1X1(t2)[1] * Der2X1(t1)[1], 2) /
                   Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 4);
        }
        else
        {
            double denominator1 = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 4);
            double numerator1 = (Der1X1(t1)[0] * (X1(t2)[0] - X1(t1)[0]) + Der1X1(t1)[1] * (X1(t2)[1] - X1(t1)[1])) *
                                (Der1X1(t2)[0] * (X1(t1)[0] - X1(t2)[0]) + Der1X1(t2)[1] * (X1(t1)[1] - X1(t2)[1]));
            double denominator2 = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
            double numerator2 = (Der1X1(t1)[0] * Der1X1(t2)[0] + Der1X1(t1)[1] * Der1X1(t2)[1]);
            
            double denominator3 = Math.Pow(Math.Sin((t1-t2)/2), 2);
            return 4 * numerator1 / denominator1 - 2 * numerator2 / denominator2 - 1 / (2 * denominator3);
        }
    }

    public static double T1J(double t, double tj, int n)
    {
        double sum = 0;
        for (int m = 1; m < n; m++)
        {
            sum += m * Math.Cos(m * (t - tj));
        }
        
        return -sum / n - Math.Cos(n * (t - tj)) / 2;
    }

    public static double T2J(double t1, double t2)
    {
        double kernelNumerator = (X2(t2)[0] - X1(t1)[0]) * VGamma1(t1)[0] + (X2(t2)[1] - X1(t1)[1]) * VGamma1(t1)[1];
        double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 2);
        
        return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
               kernelDenominator;
    }

    public static double[] GetApproximatedDerUOnGamma1(double[] uValues, int N)
    {
        double[] result = new double[2*N];
        
        for (int k = 0; k < 2*N; k++)
        {
            double tk = k * Math.PI / N;
            double s1 = 0, s2 = 0;
            for (int i = 0; i < 2*N; i++)
            {
                double ti = i * Math.PI / N;
                s2 = s2 + uValues[i + 2 * N] * T2J(tk, ti);
                s1 = s1 + (uValues[i] * T1J(tk, ti, N) + uValues[i] * HWaved(tk, ti)/(2*N)) /
                    GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0);;
            }

            result[k] = s1 + s2 / (2 * N);
        }


        return result;
    }

    public static double F1(double x)
    {
        return Math.Pow(X1(x)[0], 2) - Math.Pow(X1(x)[1], 2);
    }

    public static double G2(double x)
    {
        return 2 * X2(x)[0] * VGamma2(x)[0] - 2 * X2(x)[1] * VGamma2(x)[1];
    }

    public static double du_G1(double x)
    {
        return 2 * X1(x)[0] * VGamma1(x)[0] - 2 * X1(x)[1] * VGamma1(x)[1] ;

    }
    
    public static double GetMaxError(double[] f1Value, double[] f2Value)
    {
        double maxError = 0;

        for (int i = 0; i < f1Value.Length; i++)
        {
            maxError = Math.Max(Math.Abs(f1Value[i] - f2Value[i]), maxError);
        }

        return maxError;
    }
}

public class ND_Solver
{
    double K11(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (FunctionHelper.Der2X1(t1)[0] * FunctionHelper.Der1X1(t1)[1] - FunctionHelper.Der2X1(t1)[1] * FunctionHelper.Der1X1(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t1)[0], 0, FunctionHelper.Der1X1(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X1(t2)[0] - FunctionHelper.X1(t1)[0]) * FunctionHelper.VGamma1(t1)[0] + (FunctionHelper.X1(t2)[1] - FunctionHelper.X1(t1)[1]) * FunctionHelper.VGamma1(t1)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X1(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0) /
               kernelDenominator;
    }

    double K12(double t1, double t2)
    {
        double kernelLeftValue = (FunctionHelper.VGamma1(t1)[0] * FunctionHelper.VGamma2(t2)[0] + FunctionHelper.VGamma1(t1)[1] * FunctionHelper.VGamma2(t2)[1]) /
                                 Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X2(t2)[1]), 2);
        double kernelRightValueNominator = ((FunctionHelper.X1(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma2(t2)[0] + (FunctionHelper.X1(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma2(t2)[1]) *
                                           ((FunctionHelper.X1(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma1(t1)[0] + (FunctionHelper.X1(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma1(t1)[1]);
        double kernelRightValue = kernelRightValueNominator / 
                                  Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X2(t2)[1]), 4);
        
        return (kernelLeftValue - 2 * kernelRightValue) * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0);
    }

    double K21(double t1, double t2)
    {
        double kernelLeftValue = Math.Log(1 / FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X1(t2)[1]));
        double kernelRightValue = FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0);
        return kernelLeftValue * kernelRightValue;
    }

    double K22(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (-1 * FunctionHelper.Der1X2(t1)[0] * FunctionHelper.Der2X2(t1)[1] + FunctionHelper.Der1X2(t1)[1] * FunctionHelper.Der2X2(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t1)[0], 0, FunctionHelper.Der1X2(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X2(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma2(t2)[0] + (FunctionHelper.X2(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma2(t2)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X2(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0) /
            kernelDenominator;
    }

    public double[] solutionValues = new double[2];
    
    public double[] Solve(int N1, double[] fValues, double[] gValues)
    {
        int N = N1;
        double[] H_F_Values = new double [4*N];
        double[,] kernelMatrix = new double [4*N, 4*N];
        for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? fValues[i] : gValues[i - 2*N]; }

        for (int i = 0; i < 2*N; i++)
        {
            double ti = i * Math.PI / N;
            for (int j = 0; j < 2*N; j++)
            {
                double tj = j * Math.PI / N;
                kernelMatrix[i, j] = K11(ti, tj)/(2*N);
                kernelMatrix[i, 2*N + j] = K12(ti, tj)/(2*N);
                kernelMatrix[2*N + i, j] = K21(ti, tj)/(2*N);
                kernelMatrix[2*N + i, 2*N + j] = K22(ti, tj)/(2*N);
            }

            kernelMatrix[i, i] -= 0.5;
            kernelMatrix[2*N + i, 2*N + i] -= 0.5;
        }

        double[,] kernelMatrixExtended = new double[4*N, 4*N + 1];
        for (int i = 0; i < 4*N; i++) { kernelMatrixExtended[i, 4*N] = H_F_Values[i]; }
        for (int i = 0; i < 4*N; i++) { for (int j = 0; j < 4*N; j++) { kernelMatrixExtended[i, j] = kernelMatrix[i, j]; } }
        
        // for (int i = 0; i < 4 * N; i++) { for (int j = 0; j < 4 * N + 1; j++) { Console.Write( kernelMatrixExtended[i,j].ToString("N", setPrecision) + " "); } Console.WriteLine(); } Console.WriteLine();

        double[] ans = FunctionHelper.Gauss(kernelMatrixExtended, 4 * N);
        // Console.WriteLine("Gauss Values:");
        // for (int j = 0; j < 4 * N; j++) { Console.Write(ans[j] + " "); }
        // Console.WriteLine("\n");

        solutionValues = ans;
        return solutionValues;
    }
}

public class DN_Solver
{
    double H11(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (-1 * FunctionHelper.Der1X1(t1)[0] * FunctionHelper.Der2X1(t1)[1] + FunctionHelper.Der1X1(t1)[1] * FunctionHelper.Der2X1(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t1)[0], 0, FunctionHelper.Der1X1(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X1(t1)[0] - FunctionHelper.X1(t2)[0]) * FunctionHelper.VGamma1(t2)[0] + (FunctionHelper.X1(t1)[1] - FunctionHelper.X1(t2)[1]) * FunctionHelper.VGamma1(t2)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X1(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0) /
               kernelDenominator;
    }

    double H12(double t1, double t2)
    {
        double kernelLeftValue = Math.Log(1 / FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X2(t2)[1]));
        double kernelRightValue = FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0);
        return kernelLeftValue * kernelRightValue;
    }

    double H21(double t1, double t2)
    {
        double kernelLeftValue = (FunctionHelper.VGamma2(t1)[0] * FunctionHelper.VGamma1(t2)[0] + FunctionHelper.VGamma2(t1)[1] * FunctionHelper.VGamma1(t2)[1]) /
                                 Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X1(t2)[1]), 2);
        double kernelRightValueNominator = ((FunctionHelper.X2(t1)[0] - FunctionHelper.X1(t2)[0]) * FunctionHelper.VGamma1(t2)[0] + (FunctionHelper.X2(t1)[1] - FunctionHelper.X1(t2)[1]) * FunctionHelper.VGamma1(t2)[1]) *
                                           ((FunctionHelper.X2(t1)[0] - FunctionHelper.X1(t2)[0]) * FunctionHelper.VGamma2(t1)[0] + (FunctionHelper.X2(t1)[1] - FunctionHelper.X1(t2)[1]) * FunctionHelper.VGamma2(t1)[1]);
        double kernelRightValue = kernelRightValueNominator / 
                                  Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X1(t2)[1]), 4);
        
        return (kernelLeftValue - 2 * kernelRightValue) * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0); 
    }

    double H22(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (FunctionHelper.Der2X2(t1)[0] * FunctionHelper.Der1X2(t1)[1] - FunctionHelper.Der2X2(t1)[1] * FunctionHelper.Der1X2(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t1)[0], 0, FunctionHelper.Der1X2(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X2(t2)[0] - FunctionHelper.X2(t1)[0]) * FunctionHelper.VGamma2(t1)[0] + (FunctionHelper.X2(t2)[1] - FunctionHelper.X2(t1)[1]) * FunctionHelper.VGamma2(t1)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X2(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0) /
               kernelDenominator;
    }
    
    public double[] solutionValues = new double[2];
    
    public double[] Solve(int N1, double[] fValues, double[] gValues)
    {
        int N = N1;
        double[] H_F_Values = new double [4*N];
        double[,] kernelMatrix = new double [4*N, 4*N];
        for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? fValues[i] : gValues[i - 2*N]; }

        for (int i = 0; i < 2*N; i++)
        {
            double ti = i * Math.PI / N;
            for (int j = 0; j < 2*N; j++)
            {
                double tj = j * Math.PI / N;
                kernelMatrix[i, j] = H11(ti, tj)/(2*N);
                kernelMatrix[i, 2*N + j] = H12(ti, tj)/(2*N);
                kernelMatrix[2*N + i, j] = H21(ti, tj)/(2*N);
                kernelMatrix[2*N + i, 2*N + j] = H22(ti, tj)/(2*N);
            }
            
            kernelMatrix[i, i] += 0.5;
            kernelMatrix[2*N + i, 2*N + i] += 0.5;
        }

        double[,] kernelMatrixExtended = new double[4*N, 4*N + 1];
        for (int i = 0; i < 4*N; i++) { kernelMatrixExtended[i, 4*N] = H_F_Values[i]; }
        for (int i = 0; i < 4*N; i++) { for (int j = 0; j < 4*N; j++) { kernelMatrixExtended[i, j] = kernelMatrix[i, j]; } }
        
        // for (int i = 0; i < 4 * N; i++) { for (int j = 0; j < 4 * N + 1; j++) { Console.Write( kernelMatrixExtended[i,j].ToString("N", setPrecision) + " "); } Console.WriteLine(); } Console.WriteLine();

        double[] ans = FunctionHelper.Gauss(kernelMatrixExtended, 4 * N);
        // Console.WriteLine("Gauss Values:");
        // for (int j = 0; j < 4 * N; j++) { Console.Write(ans[j].ToString("N", setPrecision) + " "); }
        // Console.WriteLine("\n");

        
        solutionValues = ans;
        return solutionValues;
    }
}
