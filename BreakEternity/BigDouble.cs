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
    public string ToString(string? format, IFormatProvider? formatProvider) {
        throw new NotImplementedException();
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
}