package main

import (
	"fmt"
	"math"
	"math/rand"
)

xy := [][2]float64{
	{}
}

// A test eaxmple
func main() {
	xy := make([][2]float64, 10000)
	for i := 0; i < 10000; i ++ {
		r := rand.Intn(100)
		xy[i][0] = float64(i)
		xy[i][1] = math.Pow(xy[i][0], 3) - math.Pow(xy[i][0], 2) * 5 - 20 * xy[i][0] + (float64(r) - 50) / 5.0
	}
	factors := xyCurveFitting(xy, 3)
	fmt.Println(factors)
}
