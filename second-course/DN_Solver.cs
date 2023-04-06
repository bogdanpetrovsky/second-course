namespace second_course;

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
    
    public double[] Solve(int N1)
    {
        int N = N1;
        double[] H_F_Values = new double [4*N];
        double[,] kernelMatrix = new double [4*N, 4*N];
        // for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? 0 : 1; }
        for (int i = 0; i < N * 2; i++)
        {
            double ti = i * Math.PI / N;
            H_F_Values[i] = FunctionHelper.F1(ti);
            H_F_Values[i + 2*N] = FunctionHelper.G2(ti);
        }

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