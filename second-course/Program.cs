using System.Globalization;
using second_course;

NumberFormatInfo setPrecision = new NumberFormatInfo();    
setPrecision.NumberDecimalDigits = 8;
// See https://aka.ms/new-console-template for more information

void Solve(int N1)
{
    int N = N1;
    Console.WriteLine("N:" + N);
    // double[,] testCase = new double[3, 4] { { 1, 9, -5, -32 }, { -3, -5, -5, -10 }, { -2, -7, 1, 13 } }; double[] ans1 = FunctionHelper.Gauss(testCase, 3); for (int j = 0; j < 3; j++) { Console.Write(ans1[j] + " "); }
    ND_Solver ndSolver = new ND_Solver();
    DN_Solver dnSolver = new DN_Solver();
    double[] ND_Ans = ndSolver.Solve(N);
    double[] DN_Ans = dnSolver.Solve(N);
    // Console.WriteLine("ND Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(ND_Ans[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    // Console.WriteLine("DN Gauss Values:");
    // for (int j = 0; j < 4 * N; j++) { Console.Write(DN_Ans[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    // Console.WriteLine("ND ~U(1.5, 0) Value: " + FunctionHelper.GetApproximatedU(1.5, 0, ND_Ans, N, false).ToString("N", setPrecision) + " ");
    // Console.WriteLine("DN ~U(1.5, 0) Value: " + FunctionHelper.GetApproximatedU(1.5, 0, DN_Ans, N, true).ToString("N", setPrecision) + " ");

    Console.WriteLine("~U(0, 1.25) Value: " + FunctionHelper.GetApproximatedU(0, 1.25, DN_Ans, N, true).ToString("N", setPrecision) + " ");

    // double[] ND_uOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(ND_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(ND_uOnGamma1[j].ToString("N", setPrecision) + " "); }
    
    // double[] DN_uOnGamma1 = FunctionHelper.GetApproximatedUOnGamma1(DN_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(DN_uOnGamma1[j].ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
    
    double[] DN_deruOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(DN_Ans, N);
    for (int j = 0; j < 2 * N; j++) { Console.Write(Math.Abs(DN_deruOnGamma1[j]-FunctionHelper.du_G1(j*Math.PI / N)).ToString("N", setPrecision) + " "); }
    Console.WriteLine("\n");
    
    // double[] ND_deruOnGamma1 = FunctionHelper.GetApproximatedDerUOnGamma1(ND_Ans, N);
    // for (int j = 0; j < 2 * N; j++) { Console.Write(Math.Abs(ND_deruOnGamma1[j]-FunctionHelper.du_G1(j*Math.PI / N)).ToString("N", setPrecision) + " "); }
    // Console.WriteLine("\n");
}

Solve(4);
Console.WriteLine();
Solve(8);
Console.WriteLine();
Solve(16);
Console.WriteLine();
Solve(32);
Console.WriteLine();
Solve(64);
Console.WriteLine();

Console.WriteLine("\nEND");