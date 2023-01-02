using System.Text;

namespace BreakEternity;

using System;
using System.Globalization;
using Random = System.Random;

// I'm not sure if there's a "Yes, this is Unity" define symbol
// (#if UNITY doesn't seem to work). If you happen to know one - please create
// an issue here https://github.com/Razenpok/BreakInfinity.cs/issues.
#if UNITY_2017_1_OR_NEWER
using UnityEngine;
#endif
#if UNITY_2017_1_OR_NEWER
[Serializable]
#endif
public struct BigDouble : IFormattable, IComparable, IComparable<BigDouble>, IEquatable<BigDouble> {
    private const long MAX_SIGNIFICANT_DIGITS = 17; //Maximum number of digits of precision to assume in Number

    private const double
        EXP_LIMIT = 9e15; //If we're ABOVE this value, increase a layer. (9e15 is close to the largest integer that can fit in a Number.)

    private const double LAYER_DOWN = Math.Log10(9e15);

    //At layer 0, smaller non-zero numbers than this become layer 1 numbers with negative mag. After that the pattern continues as normal.
    private const double FIRST_NEG_LAYER = 1 / 9e15;

    //the largest exponent that can appear in a Double, though not all mantissas are valid here.
    private const long DOUBLE_EXP_MAX = 308;

    //The smallest exponent that can appear in a Double, though not all mantissas are valid here.
    private const long DOUBLE_EXP_MIN = -324;

    private const long MAX_ES_IN_A_ROW = 5; //For default toString behaviour, when to swap from eee... to (e^n) syntax.

    // The default size of the LRU cache used to cache Decimal.fromString.
    private const long DEFAULT_FROM_STRING_CACHE_SIZE = (1 << 10) - 1;

    private const bool IGNORE_COMMAS = true;
    private const bool COMMAS_ARE_DECIMAL_POINTS = false;

    //tetration/slog to real height stuff
//background info and tables of values for critical functions taken here: https://github.com/Patashu/break_eternity.js/issues/22
    readonly double[] critical_headers = { 2, Math.E, 3, 4, 5, 6, 7, 8, 9, 10 };

    readonly double[,] criticalTetrValues = {
        {
            // Base 2 (using http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html )
            1, 1.0891180521811202527, 1.1789767925673958433, 1.2701455431742086633, 1.3632090180450091941,
            1.4587818160364217007, 1.5575237916251418333, 1.6601571006859253673, 1.7674858188369780435,
            1.8804192098842727359,
            2,
        }, {
            // Base E (using http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html )
            1, //0.0
            1.1121114330934078681, //0.1
            1.2310389249316089299, //0.2
            1.3583836963111376089, //0.3
            1.4960519303993531879, //0.4
            1.6463542337511945810, //0.5
            1.8121385357018724464, //0.6
            1.9969713246183068478, //0.7
            2.2053895545527544330, //0.8
            2.4432574483385252544, //0.9

            Math.E, //1.0
        }, {
            // Base 3
            1, 1.1187738849693603, 1.2464963939368214, 1.38527004705667, 1.5376664685821402,
            1.7068895236551784, 1.897001227148399, 2.1132403089001035, 2.362480153784171,
            2.6539010333870774, 3,
        }, {
            // Base 4
            1, 1.1367350847096405, 1.2889510672956703, 1.4606478703324786, 1.6570295196661111,
            1.8850062585672889, 2.1539465047453485, 2.476829779693097, 2.872061932789197,
            3.3664204535587183, 4,
        }, {
            // Base 5
            1, 1.1494592900767588, 1.319708228183931, 1.5166291280087583, 1.748171114438024,
            2.0253263297298045, 2.3636668498288547, 2.7858359149579424, 3.3257226212448145,
            4.035730287722532, 5,
        }, {
            // Base 6
            1, 1.159225940787673, 1.343712473580932, 1.5611293155111927, 1.8221199554561318,
            2.14183924486326, 2.542468319282638, 3.0574682501653316, 3.7390572020926873, 4.6719550537360774,
            6,
        }, {
            // Base 7
            1, 1.1670905356972596, 1.3632807444991446, 1.5979222279405536, 1.8842640123816674,
            2.2416069644878687, 2.69893426559423, 3.3012632110403577, 4.121250340630164, 5.281493033448316,
            7,
        }, {
            // Base 8
            1, 1.1736630594087796, 1.379783782386201, 1.6292821855668218, 1.9378971836180754,
            2.3289975651071977, 2.8384347394720835, 3.5232708454565906, 4.478242031114584,
            5.868592169644505, 8,
        }, {
            // Base 9
            1, 1.1793017514670474, 1.394054150657457, 1.65664127441059, 1.985170999970283,
            2.4069682290577457, 2.9647310119960752, 3.7278665320924946, 4.814462547283592,
            6.436522247411611, 9,
        }, {
            // Base 10 (using http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html )
            1, 1.1840100246247336579, 1.4061375836156954169, 1.6802272208863963918, 2.026757028388618927,
            2.4770056063449647580, 3.0805252717554819987, 3.9191964192627283911, 5.1351528408331864230,
            6.9899611795347148455, 10,
        }
    };

    readonly double[,] critical_slog_values = {
        {
            // Base 2
            -1, -0.9194161097107025, -0.8335625019330468, -0.7425599821143978, -0.6466611521029437,
            -0.5462617907227869, -0.4419033816638769, -0.3342645487554494, -0.224140440909962,
            -0.11241087890006762, 0,
        }, {
            // Base E
            -1, //0.0
            -0.90603157029014, //0.1
            -0.80786507256596, //0.2
            -0.7064666939634, //0.3
            -0.60294836853664, //0.4
            -0.49849837513117, //0.5
            -0.39430303318768, //0.6
            -0.29147201034755, //0.7
            -0.19097820800866, //0.8
            -0.09361896280296, //0.9
            0, //1.0
        }, {
            // Base 3
            -1, -0.9021579584316141, -0.8005762598234203, -0.6964780623319391, -0.5911906810998454,
            -0.486050182576545, -0.3823089430815083, -0.28106046722897615, -0.1831906535795894,
            -0.08935809204418144, 0,
        }, {
            // Base 4
            -1, -0.8917227442365535, -0.781258746326964, -0.6705130326902455, -0.5612813129406509,
            -0.4551067709033134, -0.35319256652135966, -0.2563741554088552, -0.1651412821106526,
            -0.0796919581982668, 0,
        }, {
            // Base 5
            -1, -0.8843387974366064, -0.7678744063886243, -0.6529563724510552, -0.5415870994657841,
            -0.4352842206588936, -0.33504449124791424, -0.24138853420685147, -0.15445285440944467,
            -0.07409659641336663, 0,
        }, {
            // Base 6
            -1, -0.8786709358426346, -0.7577735191184886, -0.6399546189952064, -0.527284921869926,
            -0.4211627631006314, -0.3223479611761232, -0.23107655627789858, -0.1472057700818259,
            -0.07035171210706326, 0,
        }, {
            // Base 7
            -1, -0.8740862815291583, -0.7497032990976209, -0.6297119746181752, -0.5161838335958787,
            -0.41036238255751956, -0.31277212146489963, -0.2233976621705518, -0.1418697367979619,
            -0.06762117662323441, 0,
        }, {
            // Base 8
            -1, -0.8702632331800649, -0.7430366914122081, -0.6213373075161548, -0.5072025698095242,
            -0.40171437727184167, -0.30517930701410456, -0.21736343968190863, -0.137710238299109,
            -0.06550774483471955, 0,
        }, {
            // Base 9
            -1, -0.8670016295947213, -0.7373984232432306, -0.6143173985094293, -0.49973884395492807,
            -0.394584953527678, -0.2989649949848695, -0.21245647317021688, -0.13434688362382652,
            -0.0638072667348083, 0,
        }, {
            // Base 10
            -1, -0.8641642839543857, -0.732534623168535, -0.6083127477059322, -0.4934049257184696,
            -0.3885773075899922, -0.29376029055315767, -0.2083678561173622, -0.13155653399373268,
            -0.062401588652553186, 0,
        }
    };

    public static double decimalPlaces(double value, double places) {
        var len = places + 1;
        var numDigits = Math.Ceiling(Math.Log10(Math.Abs(value)));
        var rounded = Math.Round(value * Math.Pow(10, len - numDigits)) * Math.Pow(10, numDigits - len);
        return double.Parse(rounded.ToString("N" + Math.Max(len - numDigits, 0)));
    };

    public static double f_magLog10(double n) => Math.Sign(n) * Math.Log10(Math.Abs(n));

    //from HyperCalc source code
    public static double f_gamma(double n) {
        if (!double.IsFinite(n)) {
            return n;
        }

        if (n < -50) {
            if (n == Math.Truncate(n)) {
                return double.NegativeInfinity;
            }

            return 0;
        }

        double scal1 = 1;
        while (n < 10) {
            scal1 *= n;
            ++n;
        }

        n -= 1;
        var l = 0.9189385332046727; //0.5*Math.log(2*Math.PI)
        l += (n + 0.5) * Math.Log(n);
        l -= n;
        var n2 = n * n;
        var np = n;
        l += 1 / (12 * np);
        np *= n2;
        l += 1 / (360 * np);
        np *= n2;
        l += 1 / (1260 * np);
        np *= n2;
        l += 1 / (1680 * np);
        np *= n2;
        l += 1 / (1188 * np);
        np *= n2;
        l += 691 / (360360 * np);
        np *= n2;
        l += 7 / (1092 * np);
        np *= n2;
        l += 3617 / (122400 * np);

        return Math.Exp(l) / scal1;
    }

    const double TWO_PI = 6.2831853071795864769252842; // 2*pi
    const double EXP_N1 = 0.36787944117144232159553; // exp(-1)

    const double OMEGA = 0.56714329040978387299997; // W(1, 0)

//from https://math.stackexchange.com/a/465183
// The evaluation can become inaccurate very close to the branch point
    public double f_lambertw(double z, double tol = 1e-10) {
        double w;
        double wn;

        if (!double.IsFinite(z)) {
            return z;
        }

        switch (z) {
            case 0:
                return z;
            case 1:
                return OMEGA;
            case < 10:
                w = 0;
                break;
            default:
                w = Math.Log(z) - Math.Log(Math.Log(z));
                break;
        }

        for (var i = 0; i < 100; ++i) {
            wn = (z * Math.Exp(-w) + w * w) / (w + 1);
            if (Math.Abs(wn - w) < tol * Math.Abs(wn)) {
                return wn;
            }
            else {
                w = wn;
            }
        }

        throw new Exception($"Iteration failed to converge: ${
            z.ToString(CultureInfo.InvariantCulture)
        }");
        //return Number.NaN;
    };

//from https://github.com/scipy/scipy/blob/8dba340293fe20e62e173bdf2c10ae208286692f/scipy/special/lambertw.pxd
// The evaluation can become inaccurate very close to the branch point
// at ``-1/e``. In some corner cases, `lambertw` might currently
// fail to converge, or can end up on the wrong branch.
    public double d_lambertw(BigDouble z, double tol = 1e-10) {
        var w;
        var ew, wewz, wn;

        if (!double.IsFinite(z.mag)) {
            return z;
        }

        if (z.eq(BigDouble.dZero)) {
            return z;
        }

        if (z.eq(BigDouble.dOne)) {
            //Split out this case because the asymptotic series blows up
            return BigDouble.fromDouble(OMEGA);
        }

        //Get an initial guess for Halley's method
        w = BigDouble.ln(z);

        //Halley's method; see 5.9 in [1]

        for (let i = 0; i < 100; ++i) {
            ew = w.neg().exp();
            wewz = w.sub(z.mul(ew));
            wn = w.sub(wewz.div(w.add(1).sub(w.add(2).mul(wewz).div(Decimal.mul(2, w).add(2)))));
            if (BigDouble.abs(wn.sub(w)).lt(BigDouble.abs(wn).mul(tol))) {
                return wn;
            }
            else {
                w = wn;
            }
        }

        throw new Exception($"Iteration failed to converge: ${
            z.ToString()
        }");
        //return Decimal.dNaN;
    }

    private static BigDouble D(BigDouble value) {
        return fromValue_noAlloc(value);
    }

    private static BigDouble D(double value) {
        return fromValue_noAlloc(value);
    }

    private static BigDouble D(string value) {
        return fromValue_noAlloc(value);
    }

    private static BigDouble FC(double sign, double layer, double mag) {
        return fromComponents(sign, layer, mag);
    }

    private static BigDouble FC_NN(double sign, double layer, double mag) {
        return fromComponents_noNormalize(sign, layer, mag);
    }

    private static BigDouble ME(double mantissa, double exponent) {
        return BigDouble.fromMantissaExponent(mantissa, exponent);
    }

    private static BigDouble ME_NN(double mantissa, double exponent) {
        return BigDouble.fromMantissaExponent_noNormalize(mantissa, exponent);
    }

    public static readonly BigDouble dZero = FC_NN(0, 0, 0);
    public static readonly BigDouble dOne = FC_NN(1, 0, 1);
    public static readonly BigDouble dNegOne = FC_NN(-1, 0, 1);
    public static readonly BigDouble dTwo = FC_NN(1, 0, 2);
    public static readonly BigDouble dTen = FC_NN(1, 0, 10);
    public static readonly BigDouble dNaN = FC_NN(double.NaN, double.NaN, double.NaN);
    public static readonly BigDouble dInf = FC_NN(1, double.PositiveInfinity, double.PositiveInfinity);
    public static readonly BigDouble dNegInf = FC_NN(-1, double.NegativeInfinity, double.NegativeInfinity);
    public static readonly BigDouble dNumberMax = FC(1, 0, double.MaxValue);
    public static readonly BigDouble dNumberMin = FC(1, 0, double.MinValue);

    private static readonly LRUCache<string, BigDouble> fromStringCache = new(DEFAULT_FROM_STRING_CACHE_SIZE);

    public int sign = 0;
    public double mag = 0;
    public int layer = 0;

    public BigDouble(BigDouble value) {
        replaceFromBigDouble(value);
    }

    public BigDouble(double value) {
        replaceFromDouble(value);
    }

    public BigDouble(string value) {
        replaceFromString(value);
    }

    public double m {
        get {
            if (this.sign == 0) {
                return 0;
            }
            else if (layer == 0) {
                var exp = Math.Floor(Math.Log10(this.mag));
                //handle special case 5e-324
                double man;
                if (this.mag == 5e-324) {
                    man = 5;
                }
                else {
                    man = this.mag / powerOf10(exp);
                }

                return this.sign * man;
            }
            else if (this.layer == 1) {
                var residue = this.mag - Math.Floor(this.mag);
                return this.sign * Math.Pow(10, residue);
            }
            else {
                //mantissa stops being relevant past 1e9e15 / ee15.954
                return this.sign;
            }
        }
        set {
            if (this.layer <= 2) {
                this.fromMantissaExponent(value, this.e);
            }
            else {
                //don't even pretend mantissa is meaningful
                this.sign = Math.Sign(value);
                if (this.sign == 0) {
                    this.layer = 0;
                    this.exponent = 0;
                }
            }
        }
    }

    public double e {
        get {
            if (sign == 0) {
                return 0;
            }

            return layer switch {
                0 => Math.Floor(Math.Log10(this.mag)),
                1 => Math.Floor(this.mag),
                2 => Math.Floor(Math.Sign(this.mag) * Math.Pow(10, Math.Abs(this.mag))),
                _ => mag * double.PositiveInfinity
            };
        }
        set => this.fromMantissaExponent(this.m, value);
    }

    public double s {
        get => sign;
        set {
            if (value == 0) {
                sign = 0;
                layer = 0;
                mag = 0;
            }
            else {
                sign = value;
            }
        }
    }

    public double mantissa {
        get => m;
        set => m = value;
    }

    public double exponent {
        get => e;
        set => e = value;
    }

    public static BigDouble fromComponents(double sign, double layer, double mag) {
        return new BigDouble().replaceFromComponents(sign, layer, mag);
    }

    public static BigDouble fromComponents_noNormalize(double sign, double layer, double mag) {
        return new BigDouble().replaceFromComponentsNoNormalize(sign, layer, mag);
    }

    public static BigDouble fromMantissaExponent(double mantissa, double exponent) {
        return new BigDouble().replaceFromMantissaExponent(mantissa, exponent);
    }

    public static BigDouble fromMantissaExponent_noNormalize(double mantissa, double exponent) {
        return new BigDouble().replaceFromMantissaExponentNoNormalize(mantissa, exponent);
    }

    public static BigDouble fromBigDouble(BigDouble value) {
        return new BigDouble().replaceFromBigDouble(value);
    }

    public static BigDouble fromDouble(double value) {
        return new BigDouble().replaceFromDouble(value);
    }

    public static BigDouble fromString(string value) {
        return new BigDouble().replaceFromString(value);
    }

    public static BigDouble fromValue(BigDouble value) {
        return new BigDouble().replaceFromBigDouble(value);
    }

    public static BigDouble fromValue(double value) {
        return new BigDouble().replaceFromDouble(value);
    }

    public static BigDouble fromValue(string value) {
        return new BigDouble().replaceFromString(value);
    }

    /**
   * Converts a DecimalSource to a Decimal, without constructing a new Decimal
   * if the provided value is already a Decimal.
   *
   * As the return value could be the provided value itself, this function
   * returns a read-only Decimal to prevent accidental mutations of the value.
   * Use `new Decimal(value)` to explicitly create a writeable copy if mutation
   * is required.
   */
    private static BigDouble fromValue_noAlloc(BigDouble value) {
        return value;
    }

    private static BigDouble fromValue_noAlloc(double value) {
        return fromDouble(value);
    }

    private static BigDouble fromValue_noAlloc(string value) {
        BigDouble? cached = fromStringCache.get(value);
        if (cached != null) {
            return cached.Value;
        }

        return fromString(value);
    }

    public static BigDouble abs(dynamic value) {
        return D(value).abs();
    }

    public static BigDouble neg(dynamic value) {
        return D(value).neg();
    }

    public static BigDouble negate(dynamic value) {
        return D(value).neg();
    }

    public static BigDouble negated(dynamic value) {
        return D(value).neg();
    }

    public static BigDouble signum(dynamic value) {
        return D(value).sign;
    }

    public static BigDouble sgn(dynamic value) {
        return D(value).sign;
    }

    public static BigDouble round(dynamic value) {
        return D(value).round();
    }

    public static BigDouble floor(dynamic value) {
        return D(value).floor();
    }

    public static BigDouble ceil(dynamic value) {
        return D(value).ceil();
    }

    public static BigDouble trunc(dynamic value) {
        return D(value).trunc();
    }

    public static BigDouble add(dynamic value, dynamic other) {
        return D(value).add(other);
    }

    public static BigDouble plus(dynamic value, dynamic other) {
        return D(value).add(other);
    }

    public static BigDouble sub(dynamic value, dynamic other) {
        return D(value).sub(other);
    }

    public static BigDouble subtract(dynamic value, dynamic other) {
        return D(value).sub(other);
    }

    public static BigDouble minus(dynamic value, dynamic other) {
        return D(value).sub(other);
    }

    public static BigDouble mul(dynamic value, dynamic other) {
        return D(value).mul(other);
    }

    public static BigDouble multiply(dynamic value, dynamic other) {
        return D(value).mul(other);
    }

    public static BigDouble times(dynamic value, dynamic other) {
        return D(value).mul(other);
    }

    public static BigDouble div(dynamic value, dynamic other) {
        return D(value).div(other);
    }

    public static BigDouble divide(dynamic value, dynamic other) {
        return D(value).div(other);
    }

    public static BigDouble recip(dynamic value) {
        return D(value).recip();
    }

    public static BigDouble reciprocal(dynamic value) {
        return D(value).recip();
    }

    public static BigDouble reciprocate(dynamic value) {
        return D(value).reciprocate();
    }

    public static int cmp(dynamic value, dynamic other) {
        return D(value).cmp(other);
    }

    public static int cmpabs(dynamic value, dynamic other) {
        return D(value).cmpabs(other);
    }

    public static int compare(dynamic value, dynamic other) {
        return D(value).cmp(other);
    }

    public static bool isNaN(dynamic value) {
        value = D(value);
        return double.IsNaN(value.sign) || double.IsNaN(value.layer) || double.IsNaN(value.mag);
    }

    public static bool isFinite(dynamic value) {
        value = D(value);
        return double.IsFinite(value.sign) && double.IsFinite(value.layer) && double.IsFinite(value.mag);
    }

    public static bool eq(dynamic value, dynamic other) {
        return D(value).eq(other);
    }

    public static bool equals(dynamic value, dynamic other) {
        return D(value).eq(other);
    }

    public static bool neq(dynamic value, dynamic other) {
        return D(value).neq(other);
    }

    public static BigDouble notEquals(dynamic value, dynamic other) {
        return D(value).notEquals(other);
    }

    public static BigDouble lt(dynamic value, dynamic other) {
        return D(value).lt(other);
    }

    public static BigDouble lte(dynamic value, dynamic other) {
        return D(value).lte(other);
    }

    public static bool gt(dynamic value, dynamic other) {
        return D(value).gt(other);
    }

    public static bool gte(dynamic value, dynamic other) {
        return D(value).gte(other);
    }

    public static BigDouble max(dynamic value, dynamic other) {
        return D(value).max(other);
    }

    public static BigDouble min(dynamic value, dynamic other) {
        return D(value).min(other);
    }

    public static BigDouble minabs(dynamic value, dynamic other) {
        return D(value).minabs(other);
    }

    public static BigDouble maxabs(dynamic value, dynamic other) {
        return D(value).maxabs(other);
    }

    public static BigDouble clamp(dynamic value, dynamic min, dynamic max) {
        return D(value).clamp(min, max);
    }

    public static BigDouble clampMin(dynamic value, dynamic min) {
        return D(value).clampMin(min);
    }

    public static BigDouble clampMax(dynamic value, dynamic max) {
        return D(value).clampMax(max);
    }

    public static int cmpTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).cmp_tolerance(other, tolerance);
    }

    public static int compareTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).cmp_tolerance(other, tolerance);
    }

    public static bool eqTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).eq_tolerance(other, tolerance);
    }

    public static bool equalsTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).eq_tolerance(other, tolerance);
    }

    public static bool neqTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).neq_tolerance(other, tolerance);
    }

    public static bool notEqualsTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).notEquals_tolerance(other, tolerance);
    }

    public static bool ltTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).lt_tolerance(other, tolerance);
    }

    public static bool lteTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).lte_tolerance(other, tolerance);
    }

    public static bool gtTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).gt_tolerance(other, tolerance);
    }

    public static bool gteTolerance(dynamic value, dynamic other, double tolerance) {
        return D(value).gte_tolerance(other, tolerance);
    }

    public static BigDouble pLog10(dynamic value) {
        return D(value).pLog10();
    }

    public static BigDouble absLog10(dynamic value) {
        return D(value).absLog10();
    }

    public static BigDouble log10(dynamic value) {
        return D(value).log10();
    }

    public static BigDouble log(dynamic value, dynamic @base) {
        return D(value).log(@base);
    }

    public static BigDouble log2(dynamic value) {
        return D(value).log2();
    }

    public static BigDouble ln(dynamic value) {
        return D(value).ln();
    }

    public static BigDouble logarithm(dynamic value, dynamic @base) {
        return D(value).logarithm(@base);
    }

    public static BigDouble pow(dynamic value, dynamic other) {
        return D(value).pow(other);
    }

    public static BigDouble pow10(dynamic value) {
        return D(value).pow10();
    }

    public static BigDouble root(dynamic value, dynamic other) {
        return D(value).root(other);
    }

    public static BigDouble factorial(dynamic value) {
        return D(value).factorial();
    }

    public static BigDouble gamma(dynamic value) {
        return D(value).gamma();
    }

    public static BigDouble lngamma(dynamic value) {
        return D(value).lngamma();
    }

    public static BigDouble exp(dynamic value) {
        return D(value).exp();
    }

    public static BigDouble sqr(dynamic value) {
        return D(value).sqr();
    }

    public static BigDouble sqrt(dynamic value) {
        return D(value).sqrt();
    }

    public static BigDouble cube(dynamic value) {
        return D(value).cube();
    }

    public static BigDouble cbrt(dynamic value) {
        return D(value).cbrt();
    }

    public static BigDouble tetrate(dynamic value, double height = 2, dynamic? payload = null) {
        if (payload == null) {
            payload = FC_NN(1, 0, 1);
        }

        return D(value).tetrate(height, payload);
    }

    public static BigDouble iteratedexp(dynamic value, double height = 2, dynamic? payload = null) {
        if (payload == null) {
            payload = FC_NN(1, 0, 1);
        }

        return D(value).iteratedexp(height, payload);
    }

    public static BigDouble iteratedlog(dynamic value, dynamic? @base = null, double times = 1) {
        if (@base == null) {
            @base = 10;
        }

        return D(value).iteratedlog(@base, times);
    }

    public static BigDouble layeradd10(dynamic value, dynamic diff) {
        return D(value).layeradd10(diff);
    }

    public static BigDouble layeradd(dynamic value, double diff, double @base = 10) {
        return D(value).layeradd(diff, @base);
    }

    public static BigDouble slog(dynamic value, double @base = 10) {
        return D(value).slog(@base);
    }

    public static BigDouble lambertw(dynamic value) {
        return D(value).lambertw();
    }

    public static BigDouble ssqrt(dynamic value) {
        return D(value).ssqrt();
    }

    public static BigDouble pentate(dynamic value, double height = 2, dynamic? payload = null) {
        if (payload == null) {
            payload = FC_NN(1, 0, 1);
        }

        return D(value).pentate(height, payload);
    }

    /**
   * If you're willing to spend 'resourcesAvailable' and want to buy something
   * with exponentially increasing cost each purchase (start at priceStart,
   * multiply by priceRatio, already own currentOwned), how much of it can you buy?
   * Adapted from Trimps source code.
   */
    public static BigDouble affordGeometricSeries(
        dynamic resourcesAvailable,
        dynamic priceStart,
        dynamic priceRatio,
        dynamic currentOwned) {
        return affordGeometricSeries_core(
            D(resourcesAvailable),
            D(priceStart),
            D(priceRatio),
            currentOwned
        );
    }

    /**
   * How much resource would it cost to buy (numItems) items if you already have currentOwned,
   * the initial price is priceStart and it multiplies by priceRatio each purchase?
   */
    public static BigDouble sumGeometricSeries(
        dynamic numItems,
        dynamic priceStart,
        dynamic priceRatio,
        dynamic currentOwned) {
        return sumGeometricSeries_core(numItems, D(priceStart), D(priceRatio), currentOwned);
    }

    /**
   * If you're willing to spend 'resourcesAvailable' and want to buy something with additively
   * increasing cost each purchase (start at priceStart, add by priceAdd, already own currentOwned),
   * how much of it can you buy?
   */
    public static BigDouble affordArithmeticSeries(
        dynamic resourcesAvailable,
        dynamic priceStart,
        dynamic priceAdd,
        dynamic currentOwned) {
        return affordArithmeticSeries_core(
            D(resourcesAvailable),
            D(priceStart),
            D(priceAdd),
            D(currentOwned)
        );
    }

    /**
   * How much resource would it cost to buy (numItems) items if you already have currentOwned,
   * the initial price is priceStart and it adds priceAdd each purchase?
   * Adapted from http://www.mathwords.com/a/arithmetic_series.htm
   */
    public static BigDouble sumArithmeticSeries(
        dynamic numItems,
        dynamic priceStart,
        dynamic priceAdd,
        dynamic currentOwned) {
        return sumArithmeticSeries_core(D(numItems), D(priceStart), D(priceAdd), D(currentOwned));
    }

    /**
   * When comparing two purchases that cost (resource) and increase your resource/sec by (deltaRpS),
   * the lowest efficiency score is the better one to purchase.
   * From Frozen Cookies:
   * http://cookieclicker.wikia.com/wiki/Frozen_Cookies_(JavaScript_Add-on)#Efficiency.3F_What.27s_that.3F
   */
    public static BigDouble efficiencyOfPurchase(
        dynamic cost,
        dynamic currentRpS,
        dynamic deltaRpS) {
        return efficiencyOfPurchase_core(D(cost), D(currentRpS), D(deltaRpS));
    }

    public static BigDouble randomBigDoubleForTesting(double maxLayers) {
        var rand = new Random();
        // NOTE: This doesn't follow any kind of sane random distribution, so use this for testing purposes only.
        //5% of the time, return 0
        if (rand.NextDouble() * 20 < 1) {
            return FC_NN(0, 0, 0);
        }

        var randomsign = new Random().NextDouble() > 0.5 ? 1 : -1;

        //5% of the time, return 1 or -1
        if (rand.NextDouble() * 20 < 1) {
            return FC_NN(randomsign, 0, 1);
        }

        //pick a random layer
        var layer = Math.Floor(rand.NextDouble() * (maxLayers + 1));

        var randomexp = layer == 0 ? rand.NextDouble() * 616 - 308 : rand.NextDouble() * 16;
        //10% of the time, make it a simple power of 10
        if (rand.NextDouble() > 0.9) {
            randomexp = Math.Truncate(randomexp);
        }

        var randommag = Math.Pow(10, randomexp);
        //10% of the time, trunc mag
        if (rand.NextDouble() > 0.9) {
            randommag = Math.Truncate(randommag);
        }

        return FC(randomsign, layer, randommag);
    }

    private static BigDouble affordGeometricSeries_core(
        BigDouble resourcesAvailable,
        BigDouble priceStart,
        BigDouble priceRatio,
        dynamic currentOwned) {
        var actualStart = priceStart.mul(priceRatio.pow(currentOwned));
        return floor(
            resourcesAvailable
                .div(actualStart)
                .mul(priceRatio.sub(1))
                .add(1)
                .log10()
                .div(priceRatio.log10())
        );
    }

    private static BigDouble sumGeometricSeries_core(
        dynamic numItems,
        BigDouble priceStart,
        BigDouble priceRatio,
        dynamic currentOwned) {
        return priceStart
            .mul(priceRatio.pow(currentOwned))
            .mul(Decimal.sub(1, priceRatio.pow(numItems)))
            .div(Decimal.sub(1, priceRatio));
    }

    private static BigDouble affordArithmeticSeries_core(
        BigDouble resourcesAvailable,
        BigDouble priceStart,
        BigDouble priceAdd,
        BigDouble currentOwned) {
        // n = (-(a-d/2) + sqrt((a-d/2)^2+2dS))/d
        // where a is actualStart, d is priceAdd and S is resourcesAvailable
        // then floor it and you're done!
        var actualStart = priceStart.add(currentOwned.mul(priceAdd));
        var b = actualStart.sub(priceAdd.div(2));
        var b2 = b.pow(2);
        return b
            .neg()
            .add(b2.add(priceAdd.mul(resourcesAvailable).mul(2)).sqrt())
            .div(priceAdd)
            .floor();
    }

    private static BigDouble sumArithmeticSeries_core(
        BigDouble numItems,
        BigDouble priceStart,
        BigDouble priceAdd,
        BigDouble currentOwned) {
        var actualStart = priceStart.add(currentOwned.mul(priceAdd)); // (n/2)*(2*a+(n-1)*d)

        return numItems.div(2).mul(actualStart.mul(2).plus(numItems.sub(1).mul(priceAdd)));
    }

    public static BigDouble efficiencyOfPurchase_core(
        BigDouble cost,
        BigDouble currentRpS,
        BigDouble deltaRpS) {
        return cost.div(currentRpS).add(cost.div(deltaRpS));
    }

    public BigDouble normalize() {
        /*
        PSEUDOCODE:
        Whenever we are partially 0 (sign is 0 or mag and layer is 0), make it fully 0.
        Whenever we are at or hit layer 0, extract sign from negative mag.
        If layer === 0 and mag < FIRST_NEG_LAYER (1/9e15), shift to 'first negative layer' (add layer, log10 mag).
        While abs(mag) > EXP_LIMIT (9e15), layer += 1, mag = maglog10(mag).
        While abs(mag) < LAYER_DOWN (15.954) and layer > 0, layer -= 1, mag = pow(10, mag).
    
        When we're done, all of the following should be true OR one of the numbers is not IsFinite OR layer is not IsInteger (error state):
        Any 0 is totally zero (0, 0, 0).
        Anything layer 0 has mag 0 OR mag > 1/9e15 and < 9e15.
        Anything layer 1 or higher has abs(mag) >= 15.954 and < 9e15.
        We will assume in calculations that all Decimals are either erroneous or satisfy these criteria. (Otherwise: Garbage in, garbage out.)
        */
        if (sign == 0 || (mag == 0 && layer == 0)) {
            sign = 0;
            mag = 0;
            layer = 0;
            return this;
        }

        if (layer == 0 && this.mag < 0) {
            //extract sign from negative mag at layer 0
            mag = -mag;
            sign = -sign;
        }

        //Handle shifting from layer 0 to negative layers.
        if (layer == 0 && mag < FIRST_NEG_LAYER) {
            layer += 1;
            mag = Math.Log10(mag);
            return this;
        }

        var absmag = Math.Abs(mag);
        var signmag = Math.Sign(mag);

        if (absmag >= EXP_LIMIT) {
            layer += 1;
            mag = signmag * Math.Log10(absmag);
            return this;
        }

        while (absmag < LAYER_DOWN && layer > 0) {
            layer -= 1;
            if (layer == 0) {
                mag = Math.Pow(10, mag);
            }
            else {
                mag = signmag * Math.Pow(10, absmag);
                absmag = Math.Abs(mag);
                signmag = Math.Sign(mag);
            }
        }

        if (layer == 0) {
            if (mag < 0) {
                //extract sign from negative mag at layer 0
                mag = -mag;
                sign = -sign;
            }
            else if (this.mag == 0) {
                //excessive rounding can give us all zeroes
                sign = 0;
            }
        }

        return this;
    }

    public BigDouble replaceFromComponents(double sign, double layer, double mag) {
        this.sign = sign;
        this.layer = layer;
        this.mag = mag;

        normalize();
        return this;
    }

    public BigDouble replaceFromComponentsNoNormalize(double sign, double layer, double mag) {
        this.sign = sign;
        this.layer = layer;
        this.mag = mag;
        return this;
    }

    public BigDouble replaceFromMantissaExponent(double mantissa, double exponent) {
        layer = 1;
        sign = Math.Sign(mantissa);
        mantissa = Math.Abs(mantissa);
        mag = exponent + Math.Log10(mantissa);

        normalize();
        return this;
    }

    public BigDouble replaceFromMantissaExponentNoNormalize(double mantissa, double exponent) {
        //The idea of 'normalizing' a break_infinity.js style Decimal doesn't really apply. So just do the same thing.
        replaceFromMantissaExponent(mantissa, exponent);
        return this;
    }

    public BigDouble replaceFromBigDouble(BigDouble value) {
        sign = value.sign;
        layer = value.layer;
        mag = value.mag;
        return this;
    }

    public BigDouble replaceFromDouble(double value) {
        mag = Math.Abs(value);
        sign = Math.Sign(value);
        layer = 0;
        normalize();
        return this;
    }

    public BigDouble replaceFromString(string value) {
        var originalValue = value;
        BigDouble? cached = fromStringCache.get(originalValue);
        if (cached != null) {
            return fromBigDouble(cached.Value);
        }

        if (IGNORE_COMMAS) {
            value = value.Replace(",", "");
        }
        else if (COMMAS_ARE_DECIMAL_POINTS) {
            value = value.Replace(",", ".");
        }

        //Handle x^^^y format.
        var pentationparts = value.Split("^^^");
        double @base;
        double height;
        if (pentationparts.Length == 2) {
            @base = double.Parse(pentationparts[0]);
            height = double.Parse(pentationparts[1]);
            var heightparts = pentationparts[1].Split(";");
            var payload = 1d;
            if (heightparts.Length == 2) {
                payload = double.Parse(heightparts[1]);
                if (!isFinite(payload)) {
                    payload = 1;
                }
            }

            if (isFinite(@base) && isFinite(height)) {
                var result = pentate(base, height, payload);
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        //Handle x^^y format.
        var tetrationparts = value.Split("^^");
        if (tetrationparts.Length == 2) {
            @base = double.Parse(tetrationparts[0]);
            height = double.Parse(tetrationparts[1]);
            var heightparts = tetrationparts[1].Split(";");
            var payload = 1d;
            if (heightparts.Length == 2) {
                payload = double.Parse(heightparts[1]);
                if (!isFinite(payload)) {
                    payload = 1;
                }
            }

            if (isFinite(base) && isFinite(height)) {
                var result = tetrate(base, height, payload);
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        //Handle x^y format.
        var powparts = value.Split("^");
        double exponent;
        if (powparts.Length == 2) {
            @base = double.Parse(powparts[0]);
            exponent = double.Parse(powparts[1]);
            if (isFinite(base) && isFinite(exponent)) {
                var result = pow(base, exponent);
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        //Handle various cases involving it being a Big Number.
        value = value.Trim().ToLower(CultureInfo.InvariantCulture);

        //handle X PT Y format.
        var ptparts = value.Split("pt");
        if (ptparts.Length == 2) {
            @base = 10;
            height = double.Parse(ptparts[0]);
            ptparts[1] = ptparts[1].Replace("(", "");
            ptparts[1] = ptparts[1].Replace(")", "");
            var payload = double.Parse(ptparts[1]);
            if (!isFinite(payload)) {
                payload = 1;
            }

            if (isFinite(@base) && isFinite(height)) {
                var result = tetrate(@base, height, payload);
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        //handle XpY format (it's the same thing just with p).
        ptparts = value.Split("p");
        if (ptparts.Length == 2) {
            @base = 10;
            height = double.Parse(ptparts[0]);
            ptparts[1] = ptparts[1].Replace("(", "");
            ptparts[1] = ptparts[1].Replace(")", "");
            var payload = double.Parse(ptparts[1]);
            if (!isFinite(payload)) {
                payload = 1;
            }

            if (isFinite(@base) && isFinite(height)) {
                var result = tetrate(@base, height, payload);
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        var parts = value.Split("e");
        var ecount = parts.Length - 1;

        //Handle numbers that are exactly floats (0 or 1 es).
        if (ecount == 0) {
            var numberAttempt = double.Parse(value);
            if (isFinite(numberAttempt)) {
                replaceFromDouble(numberAttempt);
                if (fromStringCache.Count >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }
        else if (ecount == 1) {
            //Very small numbers ("2e-3000" and so on) may look like valid floats but round to 0.
            var numberAttempt = double.Parse(value);
            if (isFinite(numberAttempt) && numberAttempt != 0) {
                replaceFromDouble(numberAttempt);
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
        }

        //Handle new (e^N)X format.
        var newparts = value.Split("e^");
        if (newparts.Length == 2) {
            sign = 1;
            if (newparts[0][0] == '-') {
                sign = -1;
            }

            var layerstring = "";
            for (var i = 0; i < newparts[1].Length; ++i) {
                var chrcode = newparts[1][i];
                if ((chrcode >= 43 && chrcode <= 57) || chrcode == 101) {
                    //is "0" to "9" or "+" or "-" or "." or "e" (or "," or "/")
                    layerstring += newparts[1][i];
                } //we found the end of the layer count
                else {
                    layer = double.Parse(layerstring);
                    mag = double.Parse(newparts[1].Substring(i + 1));
                    normalize();
                    if (fromStringCache.maxSize >= 1) {
                        fromStringCache.set(originalValue, fromBigDouble(this));
                    }

                    return this;
                }
            }
        }

        if (ecount < 1) {
            sign = 0;
            layer = 0;
            mag = 0;
            if (fromStringCache.maxSize >= 1) {
                fromStringCache.set(originalValue, fromBigDouble(this));
            }

            return this;
        }

        var mantissa = double.Parse(parts[0]);
        if (mantissa == 0) {
            sign = 0;
            layer = 0;
            mag = 0;
            if (fromStringCache.maxSize >= 1) {
                fromStringCache.set(originalValue, fromBigDouble(this));
            }

            return this;
        }

        exponent = double.Parse(parts[parts.Length - 1]);
        //handle numbers like AeBeC and AeeeeBeC
        if (ecount >= 2) {
            var me = double.Parse(parts[parts.Length - 2]);
            if (isFinite(me)) {
                exponent *= Math.Sign(me);
                exponent += f_magLog10(me);
            }
        }

        //Handle numbers written like eee... (N es) X
        if (!isFinite(mantissa)) {
            sign = parts[0] == "-" ? -1 : 1;
            layer = ecount;
            mag = exponent;
        }
        //Handle numbers written like XeY
        else if (ecount == 1) {
            sign = Math.Sign(mantissa);
            layer = 1;
            //Example: 2e10 is equal to 10^log10(2e10) which is equal to 10^(10+log10(2))
            mag = exponent + Math.Log10(Math.Abs(mantissa));
        }
        //Handle numbers written like Xeee... (N es) Y
        else {
            sign = Math.Sign(mantissa);
            layer = ecount;
            if (ecount == 2) {
                var result = mul(FC(1, 2, exponent), D(mantissa));
                sign = result.sign;
                layer = result.layer;
                mag = result.mag;
                if (fromStringCache.maxSize >= 1) {
                    fromStringCache.set(originalValue, fromBigDouble(this));
                }

                return this;
            }
            else {
                //at eee and above, mantissa is too small to be recognizable!
                mag = exponent;
            }
        }

        this.normalize();
        if (fromStringCache.maxSize >= 1) {
            fromStringCache.set(originalValue, fromBigDouble(this));
        }

        return this;
    }

    public BigDouble replaceFromValue(BigDouble value) {
        return fromBigDouble(value);
    }

    public BigDouble replaceFromValue(double value) {
        return fromDouble(value);
    }

    public BigDouble replaceFromValue(string value) {
        return fromString(value);
    }


    public double toDouble() {
        if (!double.IsFinite(this.layer)) {
            return double.NaN;
        }

        return layer switch {
            0 => sign * mag,
            1 => sign * Math.Pow(10, mag),
            _ => mag > 0 ? sign > 0 ? double.PositiveInfinity : double.NegativeInfinity : 0
        };
    }

    public double mantissaWithDecimalPlaces(double places) {
        // https://stackoverflow.com/a/37425022
        if (double.IsNaN(m)) {
            return double.NaN;
        }

        if (m == 0) {
            return 0;
        }

        return decimalPlaces(this.m, places);
    }

    public double magnitudeWithDecimalPlaces(double places) {
        // https://stackoverflow.com/a/37425022
        if (double.IsNaN(mag)) {
            return double.NaN;
        }

        if (mag == 0) {
            return 0;
        }

        return decimalPlaces(this.mag, places);
    }

    public override string ToString() {
        if (double.IsNaN(layer) || double.IsNaN(sign) || double.IsNaN(mag)) {
            return "NaN";
        }

        if (double.IsPositiveInfinity(mag) || double.IsPositiveInfinity(layer)) {
            return sign == 1 ? "Infinity" : "-Infinity";
        }

        if (layer == 0) {
            if ((mag < 1e21 && mag > 1e-7) || mag == 0) {
                return (sign * mag).ToString();
            }

            return m + "e" + e;
        }

        if (layer == 1) {
            return m + "e" + e;
        }
        //layer 2+
        if (layer <= MAX_ES_IN_A_ROW) {
            return (sign == -1 ? "-" : "") + "e".repeat(layer) + mag;
        }

        return (sign == -1 ? "-" : "") + "(e^" + layer + ")" + mag;
    }

    public string toExponential(double places) {
        if (layer == 0) {
            return (sign * mag).ToString("E" + places);
        }

        return this.toStringWithDecimalPlaces(places);
    }

    public string toFixed(double places) {
        if (layer == 0) {
            return (sign * mag).ToString("N" + places);
        }

        return this.toStringWithDecimalPlaces(places);
    }

    public string toPrecision(double places) {
        if (e <= -7) {
            return toExponential(places - 1);
        }

        if (places > e) {
            return toFixed(places - exponent - 1);
        }

        return toExponential(places - 1);
    }

    public string valueOf() {
        return ToString();
    }

    public string toJSON() {
        return ToString();
    }

    public string toStringWithDecimalPlaces(double places) {
        if (layer == 0) {
            if ((mag < 1e21 && mag > 1e-7) || mag == 0) {
                return (sign * mag).ToString("N" + places);
            }

            return decimalPlaces(this.m, places) + "e" + decimalPlaces(this.e, places);
        }

        if (layer == 1) {
            return decimalPlaces(this.m, places) + "e" + decimalPlaces(this.e, places);
        }

        //layer 2+
        if (layer <= MAX_ES_IN_A_ROW) {
            return (sign == -1 ? "-" : "") + "e".repeat((int)layer) + decimalPlaces(mag, places);
        }

        return (sign == -1 ? "-" : "") + "(e^" + layer + ")" + decimalPlaces(mag, places);
    }
    
    public BigDouble abs() {
    return FC_NN(this.sign === 0 ? 0 : 1, this.layer, this.mag);
  }

  public BigDouble neg(){
    return FC_NN(-this.sign, this.layer, this.mag);
  }

  public BigDouble negate(){
    return this.neg();
  }

  public BigDouble negated(){
    return this.neg();
  }

   public signum () {
       return this.sign;
     }

  public BigDouble sgn() {
    return this.sign;
  }

  public BigDouble round() {
    if (mag < 0) {
      return dZero;
    }
    if (layer == 0) {
      return FC(sign, 0, Math.Round(mag));
    }
    return this;
  }

  public BigDouble floor() {
    if (mag < 0) {
      return dZero;
    }
    if (layer == 0) {
      return FC(sign, 0, Math.Floor(mag));
    }
    return this;
  }

  public BigDouble ceil() {
    if (mag < 0) {
      return dZero;
    }
    if (layer == 0) {
      return FC(sign, 0, Math.Ceiling(mag));
    }
    return this;
  }

  public BigDouble trunc() {
    if (mag < 0) {
      return dZero;
    }
    if (layer == 0) {
      return FC(this.sign, 0, Math.Truncate(mag));
    }
    return this;
  }

  public BigDouble add(dynamic value) {
    BigDouble dec = D(value);

    //inf/nan check
    if (!double.IsFinite(layer)) {
      return this;
    }
    if (!double.IsFinite(dec.layer)) {
      return decimal;
    }

    //Special case - if one of the numbers is 0, return the other number.
    if (sign == 0) {
      return dec;
    }
    if (dec.sign == 0) {
      return this;
    }

    //Special case - Adding a number to its negation produces 0, no matter how large.
    if (sign == -dec.sign && layer == dec.layer && mag == dec.mag) {
      return FC_NN(0, 0, 0);
    }

    BigDouble a;
    BigDouble b;

    //Special case: If one of the numbers is layer 2 or higher, just take the bigger number.
    if (layer >= 2 || dec.layer >= 2) {
      return this.maxabs(dec);
    }

    if (cmpabs(this, dec) > 0) {
      a = this;
      b = dec;
    } else {
      a = dec;
      b = this;
    }

    if (a.layer == 0 && b.layer == 0) {
      return fromDouble(a.sign * a.mag + b.sign * b.mag);
    }

    var layera = a.layer * Math.Sign(a.mag);
    var layerb = b.layer * Math.Sign(b.mag);

    //If one of the numbers is 2+ layers higher than the other, just take the bigger number.
    if (layera - layerb >= 2) {
      return a;
    }

    if (layera == 0 && layerb == -1) {
      if (Math.Abs(b.mag - Math.Log10(a.mag)) > MAX_SIGNIFICANT_DIGITS) {
        return a;
      } else {
        const magdiff = Math.pow(10, Math.log10(a.mag) - b.mag);
        const mantissa = b.sign + a.sign * magdiff;
        return FC(Math.sign(mantissa), 1, b.mag + Math.log10(Math.abs(mantissa)));
      }
    }

    if (layera === 1 && layerb === 0) {
      if (Math.abs(a.mag - Math.log10(b.mag)) > MAX_SIGNIFICANT_DIGITS) {
        return a;
      } else {
        const magdiff = Math.pow(10, a.mag - Math.log10(b.mag));
        const mantissa = b.sign + a.sign * magdiff;
        return FC(Math.sign(mantissa), 1, Math.log10(b.mag) + Math.log10(Math.abs(mantissa)));
      }
    }

    if (Math.abs(a.mag - b.mag) > MAX_SIGNIFICANT_DIGITS) {
      return a;
    } else {
      const magdiff = Math.pow(10, a.mag - b.mag);
      const mantissa = b.sign + a.sign * magdiff;
      return FC(Math.sign(mantissa), 1, b.mag + Math.log10(Math.abs(mantissa)));
    }

    throw Error("Bad arguments to add: " + this + ", " + value);
  }

  public plus(value: DecimalSource){
    return this.add(value);
  }

  public sub(value: DecimalSource){
    return this.add(D(value).neg());
  }

  public subtract(value: DecimalSource){
    return this.sub(value);
  }

  public minus(value: DecimalSource){
    return this.sub(value);
  }

  public mul(value: DecimalSource){
    const decimal = D(value);

    //inf/nan check
    if (!Number.isFinite(this.layer)) {
      return this;
    }
    if (!Number.isFinite(decimal.layer)) {
      return decimal;
    }

    //Special case - if one of the numbers is 0, return 0.
    if (this.sign === 0 || decimal.sign === 0) {
      return FC_NN(0, 0, 0);
    }

    //Special case - Multiplying a number by its own reciprocal yields +/- 1, no matter how large.
    if (this.layer === decimal.layer && this.mag === -decimal.mag) {
      return FC_NN(this.sign * decimal.sign, 0, 1);
    }

    let a;
    let b;

    //Which number is bigger in terms of its multiplicative distance from 1?
    if (
      this.layer > decimal.layer ||
      (this.layer == decimal.layer && Math.abs(this.mag) > Math.abs(decimal.mag))
    ) {
      a = this;
      b = decimal;
    } else {
      a = decimal;
      b = this;
    }

    if (a.layer === 0 && b.layer === 0) {
      return Decimal.fromNumber(a.sign * b.sign * a.mag * b.mag);
    }

    //Special case: If one of the numbers is layer 3 or higher or one of the numbers is 2+ layers bigger than the other, just take the bigger number.
    if (a.layer >= 3 || a.layer - b.layer >= 2) {
      return FC(a.sign * b.sign, a.layer, a.mag);
    }

    if (a.layer === 1 && b.layer === 0) {
      return FC(a.sign * b.sign, 1, a.mag + Math.log10(b.mag));
    }

    if (a.layer === 1 && b.layer === 1) {
      return FC(a.sign * b.sign, 1, a.mag + b.mag);
    }

    if (a.layer === 2 && b.layer === 1) {
      const newmag = FC(Math.sign(a.mag), a.layer - 1, Math.abs(a.mag)).add(
        FC(Math.sign(b.mag), b.layer - 1, Math.abs(b.mag))
      );
      return FC(a.sign * b.sign, newmag.layer + 1, newmag.sign * newmag.mag);
    }

    if (a.layer === 2 && b.layer === 2) {
      const newmag = FC(Math.sign(a.mag), a.layer - 1, Math.abs(a.mag)).add(
        FC(Math.sign(b.mag), b.layer - 1, Math.abs(b.mag))
      );
      return FC(a.sign * b.sign, newmag.layer + 1, newmag.sign * newmag.mag);
    }

    throw Error("Bad arguments to mul: " + this + ", " + value);
  }

  public multiply(value: DecimalSource){
    return this.mul(value);
  }

  public times(value: DecimalSource){
    return this.mul(value);
  }

  public div(value: DecimalSource){
    const decimal = D(value);
    return this.mul(decimal.recip());
  }

  public divide(value: DecimalSource){
    return this.div(value);
  }

  public divideBy(value: DecimalSource){
    return this.div(value);
  }

  public dividedBy(value: DecimalSource){
    return this.div(value);
  }

  public recip(){
    if (this.mag === 0) {
      return Decimal.dNaN;
    } else if (this.layer === 0) {
      return FC(this.sign, 0, 1 / this.mag);
    } else {
      return FC(this.sign, this.layer, -this.mag);
    }
  }

  public reciprocal(){
    return this.recip();
  }

  public reciprocate(){
    return this.recip();
  }

  /**
   * -1 for less than value, 0 for equals value, 1 for greater than value
   */
  public cmp(value: DecimalSource): CompareResult {
    const decimal = D(value);
    if (this.sign > decimal.sign) {
      return 1;
    }
    if (this.sign < decimal.sign) {
      return -1;
    }
    return (this.sign * this.cmpabs(value)) as CompareResult;
  }

  public cmpabs(value: DecimalSource): CompareResult {
    const decimal = D(value);
    const layera = this.mag > 0 ? this.layer : -this.layer;
    const layerb = decimal.mag > 0 ? decimal.layer : -decimal.layer;
    if (layera > layerb) {
      return 1;
    }
    if (layera < layerb) {
      return -1;
    }
    if (this.mag > decimal.mag) {
      return 1;
    }
    if (this.mag < decimal.mag) {
      return -1;
    }
    return 0;
  }

  public compare(value: DecimalSource): CompareResult {
    return this.cmp(value);
  }

  public isNan(): boolean {
    return isNaN(this.sign) || isNaN(this.layer) || isNaN(this.mag);
  }

  public isFinite(): boolean {
    return isFinite(this.sign) && isFinite(this.layer) && isFinite(this.mag);
  }

  public eq(value: DecimalSource): boolean {
    const decimal = D(value);
    return this.sign === decimal.sign && this.layer === decimal.layer && this.mag === decimal.mag;
  }

  public equals(value: DecimalSource): boolean {
    return this.eq(value);
  }

  public neq(value: DecimalSource): boolean {
    return !this.eq(value);
  }

  public notEquals(value: DecimalSource): boolean {
    return this.neq(value);
  }

  public lt(value: DecimalSource): boolean {
    return this.cmp(value) === -1;
  }

  public lte(value: DecimalSource): boolean {
    return !this.gt(value);
  }

  public gt(value: DecimalSource): boolean {
    return this.cmp(value) === 1;
  }

  public gte(value: DecimalSource): boolean {
    return !this.lt(value);
  }

  public max(value: DecimalSource){
    const decimal = D(value);
    return this.lt(decimal) ? decimal : this;
  }

  public min(value: DecimalSource){
    const decimal = D(value);
    return this.gt(decimal) ? decimal : this;
  }

  public maxabs(value: DecimalSource){
    const decimal = D(value);
    return this.cmpabs(decimal) < 0 ? decimal : this;
  }

  public minabs(value: DecimalSource){
    const decimal = D(value);
    return this.cmpabs(decimal) > 0 ? decimal : this;
  }

  public clamp(min: DecimalSource, max: DecimalSource){
    return this.max(min).min(max);
  }

  public clampMin(min: DecimalSource){
    return this.max(min);
  }

  public clampMax(max: DecimalSource){
    return this.min(max);
  }

  public cmp_tolerance(value: DecimalSource, tolerance: number): CompareResult {
    const decimal = D(value);
    return this.eq_tolerance(decimal, tolerance) ? 0 : this.cmp(decimal);
  }

  public compare_tolerance(value: DecimalSource, tolerance: number): CompareResult {
    return this.cmp_tolerance(value, tolerance);
  }

  /**
   * Tolerance is a relative tolerance, multiplied by the greater of the magnitudes of the two arguments.
   * For example, if you put in 1e-9, then any number closer to the
   * larger number than (larger number)*1e-9 will be considered equal.
   */
  public eq_tolerance(value: DecimalSource, tolerance: number): boolean {
    const decimal = D(value); // https://stackoverflow.com/a/33024979
    if (tolerance == null) {
      tolerance = 1e-7;
    }
    //Numbers that are too far away are never close.
    if (this.sign !== decimal.sign) {
      return false;
    }
    if (Math.abs(this.layer - decimal.layer) > 1) {
      return false;
    }
    // return abs(a-b) <= tolerance * max(abs(a), abs(b))
    let magA = this.mag;
    let magB = decimal.mag;
    if (this.layer > decimal.layer) {
      magB = f_maglog10(magB);
    }
    if (this.layer < decimal.layer) {
      magA = f_maglog10(magA);
    }
    return Math.abs(magA - magB) <= tolerance * Math.max(Math.abs(magA), Math.abs(magB));
  }

  public equals_tolerance(value: DecimalSource, tolerance: number): boolean {
    return this.eq_tolerance(value, tolerance);
  }

  public neq_tolerance(value: DecimalSource, tolerance: number): boolean {
    return !this.eq_tolerance(value, tolerance);
  }

  public notEquals_tolerance(value: DecimalSource, tolerance: number): boolean {
    return this.neq_tolerance(value, tolerance);
  }

  public lt_tolerance(value: DecimalSource, tolerance: number): boolean {
    const decimal = D(value);
    return !this.eq_tolerance(decimal, tolerance) && this.lt(decimal);
  }

  public lte_tolerance(value: DecimalSource, tolerance: number): boolean {
    const decimal = D(value);
    return this.eq_tolerance(decimal, tolerance) || this.lt(decimal);
  }

  public gt_tolerance(value: DecimalSource, tolerance: number): boolean {
    const decimal = D(value);
    return !this.eq_tolerance(decimal, tolerance) && this.gt(decimal);
  }

  public gte_tolerance(value: DecimalSource, tolerance: number): boolean {
    const decimal = D(value);
    return this.eq_tolerance(decimal, tolerance) || this.gt(decimal);
  }

  public pLog10(){
    if (this.lt(Decimal.dZero)) {
      return Decimal.dZero;
    }
    return this.log10();
  }

  public absLog10(){
    if (this.sign === 0) {
      return Decimal.dNaN;
    } else if (this.layer > 0) {
      return FC(Math.sign(this.mag), this.layer - 1, Math.abs(this.mag));
    } else {
      return FC(1, 0, Math.log10(this.mag));
    }
  }

  public log10(){
    if (this.sign <= 0) {
      return Decimal.dNaN;
    } else if (this.layer > 0) {
      return FC(Math.sign(this.mag), this.layer - 1, Math.abs(this.mag));
    } else {
      return FC(this.sign, 0, Math.log10(this.mag));
    }
  }

  public log(base: DecimalSource){
    base = D(base);
    if (this.sign <= 0) {
      return Decimal.dNaN;
    }
    if (base.sign <= 0) {
      return Decimal.dNaN;
    }
    if (base.sign === 1 && base.layer === 0 && base.mag === 1) {
      return Decimal.dNaN;
    } else if (this.layer === 0 && base.layer === 0) {
      return FC(this.sign, 0, Math.log(this.mag) / Math.log(base.mag));
    }

    return Decimal.div(this.log10(), base.log10());
  }

  public log2(){
    if (this.sign <= 0) {
      return Decimal.dNaN;
    } else if (this.layer === 0) {
      return FC(this.sign, 0, Math.log2(this.mag));
    } else if (this.layer === 1) {
      return FC(Math.sign(this.mag), 0, Math.abs(this.mag) * 3.321928094887362); //log2(10)
    } else if (this.layer === 2) {
      return FC(Math.sign(this.mag), 1, Math.abs(this.mag) + 0.5213902276543247); //-log10(log10(2))
    } else {
      return FC(Math.sign(this.mag), this.layer - 1, Math.abs(this.mag));
    }
  }

  public ln(){
    if (this.sign <= 0) {
      return Decimal.dNaN;
    } else if (this.layer === 0) {
      return FC(this.sign, 0, Math.log(this.mag));
    } else if (this.layer === 1) {
      return FC(Math.sign(this.mag), 0, Math.abs(this.mag) * 2.302585092994046); //ln(10)
    } else if (this.layer === 2) {
      return FC(Math.sign(this.mag), 1, Math.abs(this.mag) + 0.36221568869946325); //log10(log10(e))
    } else {
      return FC(Math.sign(this.mag), this.layer - 1, Math.abs(this.mag));
    }
  }

  public logarithm(base: DecimalSource){
    return this.log(base);
  }

  public pow(value: DecimalSource){
    const decimal = D(value);
    const a = this;
    const b = decimal;

    //special case: if a is 0, then return 0 (UNLESS b is 0, then return 1)
    if (a.sign === 0) {
      return b.eq(0) ? FC_NN(1, 0, 1) : a;
    }
    //special case: if a is 1, then return 1
    if (a.sign === 1 && a.layer === 0 && a.mag === 1) {
      return a;
    }
    //special case: if b is 0, then return 1
    if (b.sign === 0) {
      return FC_NN(1, 0, 1);
    }
    //special case: if b is 1, then return a
    if (b.sign === 1 && b.layer === 0 && b.mag === 1) {
      return a;
    }

    const result = a.absLog10().mul(b).pow10();

    if (this.sign === -1) {
      if (Math.abs(b.toNumber() % 2) % 2 === 1) {
        return result.neg();
      } else if (Math.abs(b.toNumber() % 2) % 2 === 0) {
        return result;
      }
      return Decimal.dNaN;
    }

    return result;
  }

  public pow10(){
    /*
    There are four cases we need to consider:
    1) positive sign, positive mag (e15, ee15): +1 layer (e.g. 10^15 becomes e15, 10^e15 becomes ee15)
    2) negative sign, positive mag (-e15, -ee15): +1 layer but sign and mag sign are flipped (e.g. 10^-15 becomes e-15, 10^-e15 becomes ee-15)
    3) positive sign, negative mag (e-15, ee-15): layer 0 case would have been handled in the Math.pow check, so just return 1
    4) negative sign, negative mag (-e-15, -ee-15): layer 0 case would have been handled in the Math.pow check, so just return 1
    */

    if (!Number.isFinite(this.layer) || !Number.isFinite(this.mag)) {
      return Decimal.dNaN;
    }

    let a= this;

    //handle layer 0 case - if no precision is lost just use Math.pow, else promote one layer
    if (a.layer === 0) {
      const newmag = Math.pow(10, a.sign * a.mag);
      if (Number.isFinite(newmag) && Math.abs(newmag) >= 0.1) {
        return FC(1, 0, newmag);
      } else {
        if (a.sign === 0) {
          return Decimal.dOne;
        } else {
          a = FC_NN(a.sign, a.layer + 1, Math.log10(a.mag));
        }
      }
    }

    //handle all 4 layer 1+ cases individually
    if (a.sign > 0 && a.mag >= 0) {
      return FC(a.sign, a.layer + 1, a.mag);
    }
    if (a.sign < 0 && a.mag >= 0) {
      return FC(-a.sign, a.layer + 1, -a.mag);
    }
    //both the negative mag cases are identical: one +/- rounding error
    return Decimal.dOne;
  }

  public pow_base(value: DecimalSource){
    return D(value).pow(this);
  }

  public root(value: DecimalSource){
    const decimal = D(value);
    return this.pow(decimal.recip());
  }

  public factorial(){
    if (this.mag < 0) {
      return this.add(1).gamma();
    } else if (this.layer === 0) {
      return this.add(1).gamma();
    } else if (this.layer === 1) {
      return Decimal.exp(Decimal.mul(this, Decimal.ln(this).sub(1)));
    } else {
      return Decimal.exp(this);
    }
  }

  //from HyperCalc source code
  public gamma(){
    if (this.mag < 0) {
      return this.recip();
    } else if (this.layer === 0) {
      if (this.lt(FC_NN(1, 0, 24))) {
        return Decimal.fromNumber(f_gamma(this.sign * this.mag));
      }

      const t = this.mag - 1;
      let l = 0.9189385332046727; //0.5*Math.log(2*Math.PI)
      l = l + (t + 0.5) * Math.log(t);
      l = l - t;
      const n2 = t * t;
      let np = t;
      let lm = 12 * np;
      let adj = 1 / lm;
      let l2 = l + adj;
      if (l2 === l) {
        return Decimal.exp(l);
      }

      l = l2;
      np = np * n2;
      lm = 360 * np;
      adj = 1 / lm;
      l2 = l - adj;
      if (l2 === l) {
        return Decimal.exp(l);
      }

      l = l2;
      np = np * n2;
      lm = 1260 * np;
      let lt = 1 / lm;
      l = l + lt;
      np = np * n2;
      lm = 1680 * np;
      lt = 1 / lm;
      l = l - lt;
      return Decimal.exp(l);
    } else if (this.layer === 1) {
      return Decimal.exp(Decimal.mul(this, Decimal.ln(this).sub(1)));
    } else {
      return Decimal.exp(this);
    }
  }

  public lngamma(){
    return this.gamma().ln();
  }

  public exp(){
    if (this.mag < 0) {
      return Decimal.dOne;
    }
    if (this.layer === 0 && this.mag <= 709.7) {
      return Decimal.fromNumber(Math.exp(this.sign * this.mag));
    } else if (this.layer === 0) {
      return FC(1, 1, this.sign * Math.log10(Math.E) * this.mag);
    } else if (this.layer === 1) {
      return FC(1, 2, this.sign * (Math.log10(0.4342944819032518) + this.mag));
    } else {
      return FC(1, this.layer + 1, this.sign * this.mag);
    }
  }

  public sqr(){
    return this.pow(2);
  }

  public sqrt(){
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.sqrt(this.sign * this.mag));
    } else if (this.layer === 1) {
      return FC(1, 2, Math.log10(this.mag) - 0.3010299956639812);
    } else {
      const result = Decimal.div(FC_NN(this.sign, this.layer - 1, this.mag), FC_NN(1, 0, 2));
      result.layer += 1;
      result.normalize();
      return result;
    }
  }

  public cube(){
    return this.pow(3);
  }

  public cbrt(){
    return this.pow(1 / 3);
  }

  //Tetration/tetrate: The result of exponentiating 'this' to 'this' 'height' times in a row.  https://en.wikipedia.org/wiki/Tetration
  //If payload != 1, then this is 'iterated exponentiation', the result of exping (payload) to base (this) (height) times. https://andydude.github.io/tetration/archives/tetration2/ident.html
  //Works with negative and positive real heights.
  public tetrate(height = 2, payload: DecimalSource = FC_NN(1, 0, 1)){
    //x^^1 == x
    if (height === 1) {
      return Decimal.pow(this, payload);
    }
    //x^^0 == 1
    if (height === 0) {
      return new Decimal(payload);
    }
    //1^^x == 1
    if (this.eq(Decimal.dOne)) {
      return Decimal.dOne;
    }
    //-1^^x == -1
    if (this.eq(-1)) {
      return Decimal.pow(this, payload);
    }

    if (height === Number.POSITIVE_INFINITY) {
      const this_num = this.toNumber();
      //within the convergence range?
      if (this_num <= 1.44466786100976613366 && this_num >= 0.06598803584531253708) {
        //hotfix for the very edge of the number range not being handled properly
        if (this_num > 1.444667861009099) {
          return Decimal.fromNumber(Math.E);
        }
        //Formula for infinite height power tower.
        const negln = Decimal.ln(this).neg();
        return negln.lambertw().div(negln);
      } else if (this_num > 1.44466786100976613366) {
        //explodes to infinity
        // TODO: replace this with Decimal.dInf
        return Decimal.fromNumber(Number.POSITIVE_INFINITY);
      } else {
        //0.06598803584531253708 > this_num >= 0: never converges
        //this_num < 0: quickly becomes a complex number
        return Decimal.dNaN;
      }
    }

    //0^^x oscillates if we define 0^0 == 1 (which in javascript land we do), since then 0^^1 is 0, 0^^2 is 1, 0^^3 is 0, etc. payload is ignored
    //using the linear approximation for height (TODO: don't know a better way to calculate it ATM, but it wouldn't surprise me if it's just NaN)
    if (this.eq(Decimal.dZero)) {
      let result = Math.abs((height + 1) % 2);
      if (result > 1) {
        result = 2 - result;
      }
      return Decimal.fromNumber(result);
    }

    if (height < 0) {
      return Decimal.iteratedlog(payload, this, -height);
    }

    payload = D(payload);
    const oldheight = height;
    height = Math.trunc(height);
    const fracheight = oldheight - height;

    if (this.gt(Decimal.dZero) && this.lte(1.44466786100976613366)) {
      //similar to 0^^n, flip-flops between two values, converging slowly (or if it's below 0.06598803584531253708, never. so once again, the fractional part at the end will be a linear approximation (TODO: again pending knowledge of how to approximate better, although tbh I think it should in reality just be NaN)
      height = Math.min(10000, height);
      for (let i = 0; i < height; ++i) {
        const old_payload= payload;
        payload = this.pow(payload);
        //stop early if we converge
        if (old_payload.eq(payload)) {
          return payload;
        }
      }
      if (fracheight != 0) {
        const next_payload = this.pow(payload);
        return payload.mul(1 - fracheight).add(next_payload.mul(fracheight));
      }
      return payload;
    }
    //TODO: base < 0, but it's hard for me to reason about (probably all non-integer heights are NaN automatically?)

    if (fracheight !== 0) {
      if (payload.eq(Decimal.dOne)) {
        //TODO: for bases above 10, revert to old linear approximation until I can think of something better
        if (this.gt(10)) {
          payload = this.pow(fracheight);
        } else {
          payload = Decimal.fromNumber(Decimal.tetrate_critical(this.toNumber(), fracheight));
          //TODO: until the critical section grid can handle numbers below 2, scale them to the base
          //TODO: maybe once the critical section grid has very large bases, this math can be appropriate for them too? I'll think about it
          if (this.lt(2)) {
            payload = payload.sub(1).mul(this.minus(1)).plus(1);
          }
        }
      } else {
        if (this.eq(10)) {
          payload = payload.layeradd10(fracheight);
        } else {
          payload = payload.layeradd(fracheight, this);
        }
      }
    }

    for (let i = 0; i < height; ++i) {
      payload = this.pow(payload);
      //bail if we're NaN
      if (!isFinite(payload.layer) || !isFinite(payload.mag)) {
        return payload.normalize();
      }
      //shortcut
      if (payload.layer - this.layer > 3) {
        return FC_NN(payload.sign, payload.layer + (height - i - 1), payload.mag);
      }
      //give up after 10000 iterations if nothing is happening
      if (i > 10000) {
        return payload;
      }
    }
    return payload;
  }

  //iteratedexp/iterated exponentiation: - all cases handled in tetrate, so just call it
  public iteratedexp(height = 2, payload = FC_NN(1, 0, 1)){
    return this.tetrate(height, payload);
  }

  //iterated log/repeated log: The result of applying log(base) 'times' times in a row. Approximately equal to subtracting (times) from the number's slog representation. Equivalent to tetrating to a negative height.
  //Works with negative and positive real heights.
  public iteratedlog(base: DecimalSource = 10, times = 1){
    if (times < 0) {
      return Decimal.tetrate(base, -times, this);
    }

    base = D(base);
    let result = Decimal.fromDecimal(this);
    const fulltimes = times;
    times = Math.trunc(times);
    const fraction = fulltimes - times;
    if (result.layer - base.layer > 3) {
      const layerloss = Math.min(times, result.layer - base.layer - 3);
      times -= layerloss;
      result.layer -= layerloss;
    }

    for (let i = 0; i < times; ++i) {
      result = result.log(base);
      //bail if we're NaN
      if (!isFinite(result.layer) || !isFinite(result.mag)) {
        return result.normalize();
      }
      //give up after 10000 iterations if nothing is happening
      if (i > 10000) {
        return result;
      }
    }

    //handle fractional part
    if (fraction > 0 && fraction < 1) {
      if (base.eq(10)) {
        result = result.layeradd10(-fraction);
      } else {
        result = result.layeradd(-fraction, base);
      }
    }

    return result;
  }

  //Super-logarithm, one of tetration's inverses, tells you what size power tower you'd have to tetrate base to to get number. By definition, will never be higher than 1.8e308 in break_eternity.js, since a power tower 1.8e308 numbers tall is the largest representable number.
  // https://en.wikipedia.org/wiki/Super-logarithm
  // NEW: Accept a number of iterations, and use binary search to, after making an initial guess, hone in on the true value, assuming tetration as the ground truth.
  public slog(base: DecimalSource = 10, iterations = 100){
    let step_size = 0.001;
    let has_changed_directions_once = false;
    let previously_rose = false;
    let result = this.slog_internal(base).toNumber();
    for (var i = 1; i < iterations; ++i)
    {
      let new_decimal = new Decimal(base).tetrate(result);
      let currently_rose = new_decimal.gt(this);
      if (i > 1)
      {
        if (previously_rose != currently_rose)
        {
          has_changed_directions_once = true;
        }
      }
      previously_rose = currently_rose;
      if (has_changed_directions_once)
      {
        step_size /= 2;
      }
      else
      {
        step_size *= 2;
      }
      step_size = Math.abs(step_size) * (currently_rose ? -1 : 1);
      result += step_size;
      if (step_size === 0) { break; }
    }
    return Decimal.fromNumber(result);
  }
  
  public slog_internal(base: DecimalSource = 10){
    base = D(base);

    //special cases:
    //slog base 0 or lower is NaN
    if (base.lte(Decimal.dZero)) {
      return Decimal.dNaN;
    }
    //slog base 1 is NaN
    if (base.eq(Decimal.dOne)) {
      return Decimal.dNaN;
    }
    //need to handle these small, wobbling bases specially
    if (base.lt(Decimal.dOne)) {
      if (this.eq(Decimal.dOne)) {
        return Decimal.dZero;
      }
      if (this.eq(Decimal.dZero)) {
        return Decimal.dNegOne;
      }
      //0 < this < 1: ambiguous (happens multiple times)
      //this < 0: impossible (as far as I can tell)
      //this > 1: partially complex (http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html base 0.25 for proof)
      return Decimal.dNaN;
    }
    //slog_n(0) is -1
    if (this.mag < 0 || this.eq(Decimal.dZero)) {
      return Decimal.dNegOne;
    }

    let result = 0;
    let copy = Decimal.fromDecimal(this);
    if (copy.layer - base.layer > 3) {
      const layerloss = copy.layer - base.layer - 3;
      result += layerloss;
      copy.layer -= layerloss;
    }

    for (let i = 0; i < 100; ++i) {
      if (copy.lt(Decimal.dZero)) {
        copy = Decimal.pow(base, copy);
        result -= 1;
      } else if (copy.lte(Decimal.dOne)) {
        return Decimal.fromNumber(result + Decimal.slog_critical(base.toNumber(), copy.toNumber()));
      } else {
        result += 1;
        copy = Decimal.log(copy, base);
      }
    }
    return Decimal.fromNumber(result);
  }

  //background info and tables of values for critical functions taken here: https://github.com/Patashu/break_eternity.js/issues/22
  public static slog_critical(base: number, height: number): number {
    //TODO: for bases above 10, revert to old linear approximation until I can think of something better
    if (base > 10) {
      return height - 1;
    }
    return Decimal.critical_section(base, height, critical_slog_values);
  }

  public static tetrate_critical(base: number, height: number): number {
    return Decimal.critical_section(base, height, critical_tetr_values);
  }

  public static critical_section(base: number, height: number, grid: number[][]): number {
    //this part is simple at least, since it's just 0.1 to 0.9
    height *= 10;
    if (height < 0) {
      height = 0;
    }
    if (height > 10) {
      height = 10;
    }
    //have to do this complicated song and dance since one of the critical_headers is Math.E, and in the future I'd like 1.5 as well
    if (base < 2) {
      base = 2;
    }
    if (base > 10) {
      base = 10;
    }
    let lower = 0;
    let upper = 0;
    //basically, if we're between bases, we interpolate each bases' relevant values together
    //then we interpolate based on what the fractional height is.
    //accuracy could be improved by doing a non-linear interpolation (maybe), by adding more bases and heights (definitely) but this is AFAIK the best you can get without running some pari.gp or mathematica program to calculate exact values
    //however, do note http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html can do it for arbitrary heights but not for arbitrary bases (2, e, 10 present)
    for (let i = 0; i < critical_headers.length; ++i) {
      if (critical_headers[i] == base) {
        // exact match
        lower = grid[i][Math.floor(height)];
        upper = grid[i][Math.ceil(height)];
        break;
      } else if (critical_headers[i] < base && critical_headers[i + 1] > base) {
        // interpolate between this and the next
        const basefrac =
          (base - critical_headers[i]) / (critical_headers[i + 1] - critical_headers[i]);
        lower =
          grid[i][Math.floor(height)] * (1 - basefrac) + grid[i + 1][Math.floor(height)] * basefrac;
        upper =
          grid[i][Math.ceil(height)] * (1 - basefrac) + grid[i + 1][Math.ceil(height)] * basefrac;
        break;
      }
    }
    const frac = height - Math.floor(height);
    //improvement - you get more accuracy (especially around 0.9-1.0) by doing log, then frac, then powing the result
    //(we could pre-log the lookup table, but then fractional bases would get Weird)
    //also, use old linear for slog (values 0 or less in critical section). maybe something else is better but haven't thought about what yet
    if (lower <= 0 || upper <= 0)
    {
      return lower * (1 - frac) + upper * frac;
    }
    else
    {
      return Math.pow(base, (Math.log(lower)/Math.log(base)) * (1 - frac) + (Math.log(upper)/Math.log(base)) * frac );
    }
  }

  //Function for adding/removing layers from a Decimal, even fractional layers (e.g. its slog10 representation).
  //Moved this over to use the same critical section as tetrate/slog.
  public layeradd10(diff: DecimalSource){
    diff = Decimal.fromValue_noAlloc(diff).toNumber();
    const result = Decimal.fromDecimal(this);
    if (diff >= 1) {
      //bug fix: if result is very smol (mag < 0, layer > 0) turn it into 0 first
      if (result.mag < 0 && result.layer > 0) {
        result.sign = 0;
        result.mag = 0;
        result.layer = 0;
      } else if (result.sign === -1 && result.layer == 0) {
        //bug fix - for stuff like -3.layeradd10(1) we need to move the sign to the mag
        result.sign = 1;
        result.mag = -result.mag;
      }
      const layeradd = Math.trunc(diff);
      diff -= layeradd;
      result.layer += layeradd;
    }
    if (diff <= -1) {
      const layeradd = Math.trunc(diff);
      diff -= layeradd;
      result.layer += layeradd;
      if (result.layer < 0) {
        for (let i = 0; i < 100; ++i) {
          result.layer++;
          result.mag = Math.log10(result.mag);
          if (!isFinite(result.mag)) {
            //another bugfix: if we hit -Infinity mag, then we should return negative infinity, not 0. 0.layeradd10(-1) h its this
            if (result.sign === 0) {
              result.sign = 1;
            }
            //also this, for 0.layeradd10(-2)
            if (result.layer < 0) {
              result.layer = 0;
            }
            return result.normalize();
          }
          if (result.layer >= 0) {
            break;
          }
        }
      }
    }

    while (result.layer < 0) {
      result.layer++;
      result.mag = Math.log10(result.mag);
    }
    //bugfix: before we normalize: if we started with 0, we now need to manually fix a layer ourselves!
    if (result.sign === 0) {
      result.sign = 1;
      if (result.mag === 0 && result.layer >= 1) {
        result.layer -= 1;
        result.mag = 1;
      }
    }
    result.normalize();

    //layeradd10: like adding 'diff' to the number's slog(base) representation. Very similar to tetrate base 10 and iterated log base 10. Also equivalent to adding a fractional amount to the number's layer in its break_eternity.js representation.
    if (diff !== 0) {
      return result.layeradd(diff, 10); //safe, only calls positive height 1 payload tetration, slog and log
    }

    return result;
  }

  //layeradd: like adding 'diff' to the number's slog(base) representation. Very similar to tetrate base 'base' and iterated log base 'base'.
  public layeradd(diff: number, base: DecimalSource){
    const slogthis = this.slog(base).toNumber();
    const slogdest = slogthis + diff;
    if (slogdest >= 0) {
      return Decimal.tetrate(base, slogdest);
    } else if (!Number.isFinite(slogdest)) {
      return Decimal.dNaN;
    } else if (slogdest >= -1) {
      return Decimal.log(Decimal.tetrate(base, slogdest + 1), base);
    } else {
      return Decimal.log(Decimal.log(Decimal.tetrate(base, slogdest + 2), base), base);
    }
  }

  //The Lambert W function, also called the omega function or product logarithm, is the solution W(x) === x*e^x.
  // https://en.wikipedia.org/wiki/Lambert_W_function
  //Some special values, for testing: https://en.wikipedia.org/wiki/Lambert_W_function#Special_values
  public lambertw(){
    if (this.lt(-0.3678794411710499)) {
      throw Error("lambertw is unimplemented for results less than -1, sorry!");
    } else if (this.mag < 0) {
      return Decimal.fromNumber(f_lambertw(this.toNumber()));
    } else if (this.layer === 0) {
      return Decimal.fromNumber(f_lambertw(this.sign * this.mag));
    } else if (this.layer === 1) {
      return d_lambertw(this);
    } else if (this.layer === 2) {
      return d_lambertw(this);
    }
    if (this.layer >= 3) {
      return FC_NN(this.sign, this.layer - 1, this.mag);
    }

    throw "Unhandled behavior in lambertw()";
  }

  //The super square-root function - what number, tetrated to height 2, equals this?
  //Other sroots are possible to calculate probably through guess and check methods, this one is easy though.
  // https://en.wikipedia.org/wiki/Tetration#Super-root
  public ssqrt(){
    if (this.sign == 1 && this.layer >= 3) {
      return FC_NN(this.sign, this.layer - 1, this.mag);
    }
    const lnx = this.ln();
    return lnx.div(lnx.lambertw());
  }

  //Pentation/pentate: The result of tetrating 'height' times in a row. An absurdly strong operator - Decimal.pentate(2, 4.28) and Decimal.pentate(10, 2.37) are already too huge for break_eternity.js!
  // https://en.wikipedia.org/wiki/Pentation
  public pentate(height = 2, payload: DecimalSource = FC_NN(1, 0, 1)){
    payload = D(payload);
    const oldheight = height;
    height = Math.trunc(height);
    const fracheight = oldheight - height;

    //I have no idea if this is a meaningful approximation for pentation to continuous heights, but it is monotonic and continuous.
    if (fracheight !== 0) {
      if (payload.eq(Decimal.dOne)) {
        ++height;
        payload = Decimal.fromNumber(fracheight);
      } else {
        if (this.eq(10)) {
          payload = payload.layeradd10(fracheight);
        } else {
          payload = payload.layeradd(fracheight, this);
        }
      }
    }

    for (let i = 0; i < height; ++i) {
      payload = this.tetrate(payload.toNumber());
      //bail if we're NaN
      if (!isFinite(payload.layer) || !isFinite(payload.mag)) {
        return payload.normalize();
      }
      //give up after 10 iterations if nothing is happening
      if (i > 10) {
        return payload;
      }
    }

    return payload;
  }

  // trig functions!
  public sin(): this | Decimal {
    if (this.mag < 0) {
      return this;
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.sin(this.sign * this.mag));
    }
    return FC_NN(0, 0, 0);
  }

  public cos(){
    if (this.mag < 0) {
      return Decimal.dOne;
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.cos(this.sign * this.mag));
    }
    return FC_NN(0, 0, 0);
  }

  public tan(): this | Decimal {
    if (this.mag < 0) {
      return this;
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.tan(this.sign * this.mag));
    }
    return FC_NN(0, 0, 0);
  }

  public asin(): this | Decimal {
    if (this.mag < 0) {
      return this;
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.asin(this.sign * this.mag));
    }
    return FC_NN(Number.NaN, Number.NaN, Number.NaN);
  }

  public acos(){
    if (this.mag < 0) {
      return Decimal.fromNumber(Math.acos(this.toNumber()));
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.acos(this.sign * this.mag));
    }
    return FC_NN(Number.NaN, Number.NaN, Number.NaN);
  }

  public atan(): this | Decimal {
    if (this.mag < 0) {
      return this;
    }
    if (this.layer === 0) {
      return Decimal.fromNumber(Math.atan(this.sign * this.mag));
    }
    return Decimal.fromNumber(Math.atan(this.sign * 1.8e308));
  }

  public sinh(){
    return this.exp().sub(this.negate().exp()).div(2);
  }

  public cosh(){
    return this.exp().add(this.negate().exp()).div(2);
  }

  public tanh(){
    return this.sinh().div(this.cosh());
  }

  public asinh(){
    return Decimal.ln(this.add(this.sqr().add(1).sqrt()));
  }

  public acosh(){
    return Decimal.ln(this.add(this.sqr().sub(1).sqrt()));
  }

  public atanh(){
    if (this.abs().gte(1)) {
      return FC_NN(Number.NaN, Number.NaN, Number.NaN);
    }

    return Decimal.ln(this.add(1).div(Decimal.fromNumber(1).sub(this))).div(2);
  }

  /**
   * Joke function from Realm Grinder
   */
  public ascensionPenalty(ascensions: DecimalSource){
    if (ascensions === 0) {
      return this;
    }

    return this.root(Decimal.pow(10, ascensions));
  }

  /**
   * Joke function from Cookie Clicker. It's 'egg'
   */
  public egg(){
    return this.add(9);
  }

  public lessThanOrEqualTo(other: DecimalSource): boolean {
    return this.cmp(other) < 1;
  }

  public lessThan(other: DecimalSource): boolean {
    return this.cmp(other) < 0;
  }

  public greaterThanOrEqualTo(other: DecimalSource): boolean {
    return this.cmp(other) > -1;
  }

  public greaterThan(other: DecimalSource): boolean {
    return this.cmp(other) > 0;
  }

  public int CompareTo(object? obj) {
        throw new NotImplementedException();
    }

    public int CompareTo(BigDouble other) {
        throw new NotImplementedException();
    }

    public bool Equals(BigDouble other) {
        throw new NotImplementedException();
    }

    /// <summary>
    /// We need this lookup table because Math.pow(10, exponent) when exponent's absolute value
    /// is large is slightly inaccurate. you can fix it with the power of math... or just make
    /// a lookup table. Faster AND simpler.
    /// </summary>
    private static class PowersOf10 {
        private static double[] POWERS { get; } = new double[DOUBLE_EXP_MAX - DOUBLE_EXP_MIN];
        private const long INDEX_OF_0 = -DOUBLE_EXP_MIN - 1;

        static PowersOf10() {
            var index = 0;
            for (var i = 0; i < POWERS.Length; i++) {
                POWERS[index++] = double.Parse("1e" + (i - INDEX_OF_0), CultureInfo.InvariantCulture);
            }
        }

        public static double lookup(long power) {
            return POWERS[INDEX_OF_0 + power];
        }
    }

    private struct PrivateConstructorArg {
    }
}

public static class StringExtensions {
    public static string repeat(this string s, int n)
        => new StringBuilder(s.Length * n).Insert(0, s, n).ToString();
}