using System.Globalization;
using second_course;

NumberFormatInfo setPrecision = new NumberFormatInfo();    
setPrecision.NumberDecimalDigits = 8;

void Solve(int N1)
{
    int N = N1;
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

        g1Values[i] = 0;
        f2Values[i] = FunctionHelper.F2(ti, false);
        g2Values[i] = FunctionHelper.G2(ti, false);
    }
    
    // double[,] testCase = new double[3, 4] { { 1, 9, -5, -32 }, { -3, -5, -5, -10 }, { -2, -7, 1, 13 } }; double[] ans1 = FunctionHelper.Gauss(testCase, 3); for (int j = 0; j < 3; j++) { Console.Write(ans1[j] + " "); }
    ND_Solver ndSolver = new ND_Solver();
    DN_Solver dnSolver = new DN_Solver();
    
    double solutionError = Double.MaxValue;

    for (int ii = 0; ii < 20; ii++)
    {
        double[] ND_Ans = ndSolver.Solve(N, g1Values, f2Values);
        double[] ND_UOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(ND_Ans, N); 
        if (ii == 0)
        { 
            f1Values = ND_UOnGamma1;
        }
        else
        {
            
            for (int jj = 0; jj < 2 * N; jj++)
            {
                f1Values[jj] = FunctionHelper.relexationParameter * ND_UOnGamma1[jj] + (1 - FunctionHelper.relexationParameter) * f1Values[jj];
            }
        }
        
        
        // Console.WriteLine("F1 Values:");
        // for (int j = 0; j < f1Values.Length; j++)
        // {
            // Console.WriteLine(f1Values[j] + " ");
        // }
        
        // Console.WriteLine("G1 Values:");
        // for (int j = 0; j < g2Values.Length; j++)
        // {
        //     Console.WriteLine(g2Values[j] + " ");
        // }

        double[] DN_Ans = dnSolver.Solve(N, f1Values, g2Values);
        double[] DN_derUOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(DN_Ans, N);
        g1Values = DN_derUOnGamma1;

        // Console.WriteLine("G1 Values:");
        // for (int j = 0; j < g1Values.Length; j++)
        // {
            // Console.WriteLine(g1Values[j] + " ");
        // }
        
        Console.WriteLine("Solution Error:");
        solutionError = FunctionHelper.GetMaxError(f1Values, exactSolution);
        // double maxX = -1;
        // for (int i = 0; i < N*2; i++)
        // {
            // maxX = Math.Max(g1Values[i], maxX);
        // }
        
        Console.WriteLine(solutionError);
        // Console.WriteLine(maxX);
    }
    
    // while (solutionError > FunctionHelper.ErrorEps)
    // {
    //     double[] DN_Ans = dnSolver.Solve(N, f1Values, g2Values);
    //     double[] DN_derUOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(DN_Ans, N);
    //     g1Values = DN_derUOnGamma1;
    //     
    //     double[] ND_Ans = ndSolver.Solve(N, f2Values, g1Values);
    //     double[] ND_UOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(ND_Ans, N);
    //     f1Values = ND_UOnGamma1;
    //
    //     solutionError = FunctionHelper.GetMaxError(f1Values, exactSolution);
    //     Console.WriteLine(solutionError);
    // }
    
    
    
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

// Solve(4);
// Console.WriteLine();
// Solve(8);
// Console.WriteLine();
// Solve(16);
// Console.WriteLine();
// Solve(32);
// Console.WriteLine();
Solve(64);
Console.WriteLine();

Console.WriteLine("\nEND");

