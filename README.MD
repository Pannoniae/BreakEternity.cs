# BreakEternity.cs
A low-quality C# port of [break_eternity.js](https://github.com/Patashu/break_eternity.js) - a solution for incremental gameswhich want to represent numbers as large as 10^^1e308 and as 'small' as 10^-(10^^1e308). This is a sequel to break_infinity.js, the other library which maxes out at 1e1e308 ( https://github.com/Patashu/break_infinity.js ) and its C# port ( https://github.com/Razenpok/BreakInfinity.cs ).  
Despite handling a wider range of numbers, execution time is comparable (within 2x/0.5x as fast as break_infinity.js in testing) and it has the same interface, so it can be used as a drop-in replacement for break_infinity.js and decimal.js.which want to deal with very large numbers (bigger in magnitude than 1e308, up to as much as 1e(9e15) ) and want to prioritize speed over accuracy.

## Installation

Just drop [BigDouble.cs](https://github.com/Pannoniae/BreakEternity.cs/blob/master/BreakEternity/BigDouble.cs) file into your `Scripts` folder and use `BigDouble` type instead of `double` in your scripts

## BigDouble
`BigDouble` is a `double` replacement for very large numbers.

# Special thanks:
- [Patashu](https://github.com/Patashu) - for amazing library and major effort in porting to C#
- https://mrob.com/pub/comp/hypercalc/hypercalc-javascript.html HyperCalc, an existing calculator that handles numbers until 10^^(1e10) and is proving very useful for testing. (I also use SpeedCrunch, which goes up to 1e1e9, and break_infinity.js, which goes up to 1e1e308, for testing)
- https://play.google.com/store/apps/details?id=com.antoine.mathematician.oddlittlegame&hl=en Incremental Unlimited, an incremental game that reaches as high as 10^^4.
- http://myweb.astate.edu/wpaulsen/tetcalc/tetcalc.html Tetration Calculator, a javascript imeplementation of tetration to continuous heights (with pre-computed bases) with high precision.
- All the folks over at Tetration Forum working slowly but surely to make calculators for our favourite non-elementary functions: https://math.eretrandre.org/tetrationforum/forumdisplay.php?fid=8
- nathanisbored, for coming up with a short and sweet formula for 10^a + 10^b == 10^c which has been used in add and mul.
- Razenpok, slabdrill and Hevipelle/Antimatter Dimensions, for inspiration, assistance and testing with the original break_infinity.js and its C# port.

obligatory SEO: number library, big number, big num, bignumber, bignum, big integer, biginteger, bigint, incremental games, idle games, large numbers, huge numbers

# Talk to me, I'm lonely

https://discord.gg/wcqwNPpFqP
