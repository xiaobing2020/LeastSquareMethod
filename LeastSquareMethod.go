package main

import "math"

func isZero(zeros []float64) bool {
	const zero = 1e-12
	for i := 0; i < len(zeros); i++ {
		if zeros[i] > zero {
			return false
		}
	}
	return true
}

// 计算矩阵第index行 被 matrix[i][i]分割的两部分的和(带factors系数)
func calculateSeparateByIndex(matrix [][]float64,factors []float64,
	order int, index int) (lessIndexPart,biggerIndexPart float64) {
	lessIndexPart = 0
	biggerIndexPart = 0
    for i := 0; i < index; i++ {
		lessIndexPart += factors[i] * matrix[index][i]
	}
	for i := index + 1; i < order + 1; i++ {
		biggerIndexPart += factors[i] * matrix[index][i]
	}
    return lessIndexPart, biggerIndexPart
}

func calculateFactors(matrix [][]float64, order int) (factors []float64) {
	factors = make([]float64, order + 1)
	factorsTemp := make([]float64, order + 1)
	factorsZero := make([]float64, order + 1)
    for i := 0; i < order + 1; i++ {
		factorsZero[i] = 1.0
	}

	// Jacobi迭代计算拟合函数的系数
	for !isZero(factorsZero) {
		for i := 0; i < order + 1; i++ {
			factorsTemp[i] = factors[i]
			A,B := calculateSeparateByIndex(matrix, factors, order, i)
			factors[i] = (matrix[i][order + 1] - A - B) / matrix[i][i]

			// 迭代条件更新
			factorsZero[i] = (factors[i] - factorsTemp[i]) * (factors[i] - factorsTemp[i])
		}
	}
	return factors
}

// xyCurveFitting: 最小二乘法计算拟合曲线
//       参数  xyArray: 原始xy数据对
//       参数    order: 预期拟合的曲线最高阶次数
//       返回值 factors: 拟合曲线各项系数 从常数项系数 ~ order次项系数
func xyCurveFitting(xyArray [][2]float64, order int) (factors []float64) {
	dataLength := len(xyArray)

	// 计算 ∑xi， ∑xi^2，... ~... ∑xi^2n
	xin := make([]float64, order * 2 + 1) // ∑(xi^n)
	xin[0] = float64(dataLength)
	for i := 1; i < order * 2 + 1; i++ {
		for j := 0; j < dataLength; j++ {
			xin[i] += math.Pow(xyArray[j][0], float64(i))
		}
	}

	// 计算 ∑yi, ∑xiyi, ∑xi^2*yi,... ~,... ∑ xi^n * yi
    xinyi := make([]float64, order + 1)  // ∑(xi^n * yi)
	for i := 0; i < order + 1; i++ {
		for j := 0; j < dataLength; j++ {
			if i == 0 {
				xinyi[i] += xyArray[j][1]
			} else {
				xinyi[i] += math.Pow(xyArray[j][0], float64(i)) * xyArray[j][1]
			}
		}
	}

	// 创建 (order+1)*(order+2)增广矩阵存放待解方程组系数
	matrix := make([][]float64, order + 1)
    for  i := 0; i < order + 1; i++ {
 	   matrix[i] = make([]float64, order + 2)
    }


    //  | dataLenth ∑xi ∑xi^2 ...      ∑xi^n    ∑yi   |
    //  | ∑xi ∑xi^2 ∑xi^3 ...      ∑xi^n+1  ∑xi*yi |
    //  | ... ... ... ... .... ... .. ..... ... ..|
    //  | ∑xi^n ∑xi^n+1 ∑xi^n+2 ...∑xi^2n   ∑xi^n *yi |
	// 增广矩阵初始化
	for i := 0; i < order + 1; i++ {
		matrix[i][order + 1] = xinyi[i]
	}
	for i := 0; i < order + 1; i++ {
		for j := 0; j < order + 1; j++ {
			matrix[i][j] = xin[i + j]
		}
	}
	
	// 迭代计算拟合曲线系数
	factors = calculateFactors(matrix, order)

	return factors
}
