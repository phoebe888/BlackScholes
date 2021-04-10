/******************************************************************************

                            Online C# Compiler.
                Code, Compile, Run and Debug C# program online.
Write your code in this editor and press "Run" button to execute it.

*******************************************************************************/
using System;
class Program
{
    static void Main()
    {
        double K;
        double S;
        double r;
        double q;
        double vol;
        double t;
        double t0;
        string optype;
        double target;
        double impvol = 0.01;
        Console.WriteLine("please enter an option strike price = ");
        K = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter a current Stock price = ");
        S = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter a risk free rate(%) = ");
        r = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter a dividend yield(%) = ");
        q = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter a stock volatility(%) = ");
        vol = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter T, expiration time = ");
        t = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter t0, start time = ");
        t0 = double.Parse(Console.ReadLine());
        Console.WriteLine("please enter the option: call or put");
        optype = Console.ReadLine();
        RandLib rd = new RandLib();
        BSOption option = new BSOption(K, r, q, vol, t, optype);
        double fairvalue = BSFairValue(S, t0, option);
        double delta = Delta(S, t0, option);
        double gamma = Gamma(S, t0, option);
        double vega = Vega(S, t0, option);
        double rho = Rho(S, t0, option);
        double theta = Theta(S, t0, option);
        Console.WriteLine($"The fair value of the option is {fairvalue}");
        Console.WriteLine($"The delta of the option is {delta}");
        Console.WriteLine($"The gamma of the option is {gamma}");
        Console.WriteLine($"The vega of the option is {vega}");
        Console.WriteLine($"The rho of the option is {rho}");
        Console.WriteLine($"The theta of the option is {theta}");
        Console.WriteLine("please enter the target option value to calculate implied volatility");
        target = double.Parse(Console.ReadLine());
        BSImpvol(S, t0, target, 5, 0.001, 0.001, ref impvol, option);
        Console.WriteLine($"The implied volatility of the option is {impvol}");
    }
    
    public static double BSFairValue(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double value;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            value = 0;
        }
        if (opttype == "call")
        {
            value = Math.Round((S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(d1) - K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(d2)), 4);
        }
        else
        {
            value = Math.Round((K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(-d2) - S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(-d1)), 4);
        }
        return value;
    }
    
    public static double Delta(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double dvalue;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            dvalue = 0;
        }
        if (opttype == "call")
        {
            dvalue = Math.Round(Math.Exp(-q * (T - t0)) * (RandLib.NormalCDF((Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0)))), 4);
        }
        else
        {
            dvalue = Math.Round(Math.Exp(-q * (T - t0)) * (RandLib.NormalCDF((Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0)))) - 1, 4);
        }
        return dvalue;
    }
    
    public static double Gamma(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double gvalue;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            gvalue = 0;
        }
        gvalue = Math.Round((Math.Exp(-q * (T - t0)) * RandLib.NormalPDF(d1)) / (S * vol * Math.Sqrt(T - t0)), 4);
        return gvalue;
    }
    
    public static double Vega(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double vvalue;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            vvalue = 0;
        }
        vvalue = Math.Round(.01 * S * Math.Exp(-q * (T - t0)) * Math.Sqrt(T - t0) * RandLib.NormalPDF(d1), 4);
        return vvalue;
    }
    
    public static double Rho(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double rvalue;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            rvalue = 0;
        }
        if (opttype == "call")
        {
            rvalue = Math.Round(.01 * (T - t0) * K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(d2), 4);
        }
        else
        {
            rvalue = Math.Round(-.01 * (T - t0) * K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(-d2), 4);
        }
        return rvalue;
    }
    
    public static double Theta(double S, double t0, BSOption opt)
    {
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        string opttype = opt.getoption();
        double thetavalue;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);
        if (S <= 0 || K <= 0 || vol <= 0 || T - t0 <= 0)
        {
            thetavalue = 0;
        }
        if (opttype == "call")
        {
            thetavalue = Math.Round(((-S * Math.Exp(-q * (T - t0)) * RandLib.NormalPDF(d1) * 0.5 * vol / (Math.Sqrt(T - t0))) + q * S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(d1) - r * K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(d2)) / 365, 4);
        }
        else
        {
            thetavalue = Math.Round((-S * Math.Exp(-q * (T - t0)) * RandLib.NormalPDF(d1) * 0.5 * vol / (Math.Sqrt(T - t0)) + r * K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(-d2) - q * S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(-d1)) / 365, 4);
        }
        return thetavalue;
    }
    
    public static bool BSImpvol(double S, double t0, double target, int max_iter, double Opt_tol, double vol_tol, ref double impvol, BSOption opt)
    {
        bool bool_value = true;
        double T = opt.gettimeexp();
        double K = opt.getstrike();
        double vol = opt.getvolatility();
        double q = opt.getyield();
        double r = opt.getriskfreerate();
        double vol_in = 0.4;
        double vol_f = 0;
        string opttype = opt.getoption();
        double V;
        double vega;
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * Math.Pow(vol, 2)) * (T - t0)) / (vol * Math.Sqrt(T - t0));
        double d2 = d1 - vol * Math.Sqrt(T - t0);

        for (int i = 0; i <= max_iter; i++)
        {
            if (opttype == "call")
            {
                V = S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(d1) - K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(d2);
            }
            else
            {
                V = K * Math.Exp(-r * (T - t0)) * RandLib.NormalCDF(-d2) - S * Math.Exp(-q * (T - t0)) * RandLib.NormalCDF(-d1);
            }
            vega = S * Math.Exp(-q * (T - t0)) * RandLib.NormalPDF(d1) * Math.Sqrt(T - t0);
            if (i == 0)
            {
                vol_f = vol_in - (V - target) / vega;
            }
            else
            {
                if (Math.Abs(V - target) <= Opt_tol && Math.Abs((vol_f - vol_in)) <= vol_tol)
                {
                    vol_f = vol_in;
                    vol_f = vol_in - (V - target) / vega;
                }
            }
        }
        impvol = vol_f;
        return bool_value;
    }
}

class BSOption
{
    private double k;
    private double r;
    private double q;
    private double T;
    private string putcall;
    private double vol;
    public BSOption(double strike, double rfrate, double yield, double volatility, double time, string option)
    {
        k = strike;
        r = rfrate / 100;
        q = yield / 100;
        T = time;
        putcall = option;
        vol = volatility / 100;
    }
    public double getstrike() { return k; }
    public double getriskfreerate() { return r; }
    public double getyield() { return q; }
    public double gettimeexp() { return T; }
    public double getvolatility() { return vol; }
    public string getoption() { return putcall; }
    public void print()
    {
        Console.WriteLine($"Strike = {k}");
        Console.WriteLine($"Risk Free Rate = {r}");
        Console.WriteLine($"Yield = {q}");
        Console.WriteLine($"Volatility = {vol}");
        Console.WriteLine($"Expiration Time = {T}");
        Console.WriteLine($"The option is {putcall}");
    }
}


sealed class RandLib
{
    private Random r;
    public RandLib()
    {
        r = new Random();
    }
    public double NextDouble()
    {
        return r.NextDouble();
    }
    public double NextDouble(double b)
    {
        if (b == 0.0) return 0.0;
        if (b < 0.0) return NextDouble(b, 0.0);
        return b * r.NextDouble();
    }
    public double NextDouble(double a, double b)
    {
        if (b == a) return a;
        if (b < a) return NextDouble(b, a);
        return a + (b - a) * r.NextDouble();
    }
    public double NextExponential()
    {
        double x = 1.0 - r.NextDouble();
        return - Math.Log(x);
    }
    public double NextNormal()
    {
        double U = 1.0;
        double V = 1.0;
        double S = (U * U) + (V * V);
        while ((S > 1.0) || (S == 0.0))
        {
            U = 2.0 * r.NextDouble() - 1.0;
            V = 2.0 * r.NextDouble() - 1.0;
            S = (U * U) + (V * V);
        }
        double tmp = Math.Sqrt(-2.0 * Math.Log(S) / S);
        double X = U * tmp;
        double Y = V * tmp;
        return X;
    }
    public double NextNormal(double mu, double sigma)
    {
        return (mu + sigma * NextNormal());
    }
    public static double ExponentialPDF(double x)
    {
        if (x <= 0.0) return 0.0;
        return Math.Exp(-x);
    }
    public static double ExponentialCDF(double x)
    {
        if (x <= 0.0) return 0.0;
        return 1.0 - Math.Exp(-x);
    }
    public static double NormalPDF(double x)
    {
        return Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);
    }
    public static double NormalCDF(double x)
    {
        if (x == 0.0) return 0.5;
        if (x < 0.0) return 1.0 - NormalCDF(-x);
        double[] b = { 0.2316419, 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429 };
        double t = 1.0 / (1.0 + b[0] * x);
        double sum = t * (b[1] + t * (b[2] + t * (b[3] + t * (b[4] + t * (b[5])))));
        double Phi = 1.0 - NormalPDF(x) * sum;
        return Phi;
    }
}

