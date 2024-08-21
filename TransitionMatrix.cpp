#include "TransitionMatrix.hpp"
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <format>
#include <bit>

TransitionMatrix::TransitionMatrix()
{
}

TransitionMatrix::~TransitionMatrix()
{
}

void TransitionMatrix::SetIdentity(uint32_t gridWidth, uint32_t gridHeight)
{
	uint8_t identityRule = 0b010;
	SetRuleSquareGrid(gridWidth, gridHeight, boost::dynamic_bitset<uint64_t>(3, identityRule), 3, 1);
}

TransitionMatrix TransitionMatrix::Mul(const TransitionMatrix& right)
{
	if(right.mGridWidth != mGridWidth || right.mGridHeight != mGridHeight)
	{
		return *this;
	}

	TransitionMatrix result;
	result.mGridWidth  = mGridWidth;
	result.mGridHeight = mGridHeight;

	result.mRows.resize(mRows.size());
	for(size_t i = 0; i < mRows.size(); i++)
	{
		result.mRows[i].resize(mRows[i].size(), 0);
	}

	assert(right.mRows.size() > 0 && right.mRows.size() == right.mRows[0].size());
	
	std::vector<boost::dynamic_bitset<uint64_t>> rightColumnMajor;
	rightColumnMajor.resize(right.mRows.size());
	for(size_t j = 0; j < right.mRows[0].size(); j++)
	{
		rightColumnMajor[j].resize(right.mRows.size());
	}

	for(size_t i = 0; i < right.mRows.size(); i++)
	{
		for(size_t j = 0; j < mRows.size(); j++)
		{
			rightColumnMajor[j].set(i, right.mRows[i][j]);
		}
	}

#pragma omp parallel for
	for(int i = 0; i < (int)mRows.size(); i++)
	{
		const boost::dynamic_bitset<uint64_t>& leftMatrixRow = mRows[i];

		for(size_t j = 0; j < result.mRows[i].size(); j++)
		{
			const boost::dynamic_bitset<uint64_t>& rightMatrixColumn = rightColumnMajor[j];

			bool resMul = ((leftMatrixRow & rightMatrixColumn).count() % 2) == 1;
			result.mRows[i][j] = resMul;
		}
	}

	return result;
}

boost::dynamic_bitset<uint64_t> TransitionMatrix::MulState(const StateVector& state)
{
	assert(state.size() == mGridWidth * mGridHeight);
	boost::dynamic_bitset<uint64_t> result(state.size());

	for(uint32_t i = 0; i < mGridWidth * mGridHeight; i++)
	{
		result.set(i, (mRows[i] & state).count() % 2);
	}

	return result;
}

TransitionMatrix TransitionMatrix::CalcMatrixPower(const boost::multiprecision::cpp_int& matrixPower)
{
	TransitionMatrix mulMat;

	if(matrixPower == 1)
	{
		mulMat = *this;
	}
	else
	{
		boost::multiprecision::cpp_int totalMatrixPower = 0;
		mulMat.SetIdentity(mGridWidth, mGridHeight);

		if(matrixPower & 1)
		{
			mulMat = mulMat.Mul(*this);
		}

		boost::multiprecision::cpp_int currMatrixPower = 1;

		//For each "1" bit of matrixPower, multiply the result matrix by the current matrix
		TransitionMatrix prevMatrix = *this;
		boost::multiprecision::cpp_int currMatrixPowerBit = 2; //Power 1 is already checked
		while(currMatrixPowerBit <= matrixPower)
		{
			TransitionMatrix currMatrix = prevMatrix.Mul(prevMatrix);
			if(matrixPower & currMatrixPowerBit)
			{
				mulMat = mulMat.Mul(currMatrix);

				currMatrixPower = currMatrixPower | currMatrixPowerBit;
			}

			std::swap(prevMatrix, currMatrix);
			currMatrixPowerBit = (currMatrixPowerBit << 1);
		}
	}

	return mulMat;
}

void TransitionMatrix::SetRule150Square(uint32_t gridWidth)
{
	uint8_t rule150Bits = 0b111;
	SetRuleSquareGrid(gridWidth, 1, boost::dynamic_bitset<uint64_t>(3, rule150Bits), 3, 1);
}

void TransitionMatrix::SetRule150Cyclic(uint32_t gridWidth)
{
	uint8_t rule150Bits = 0b111;
	SetRuleToroidGrid(gridWidth, 1, boost::dynamic_bitset<uint64_t>(3, rule150Bits), 3, 1);
}

void TransitionMatrix::SetRule90Square(uint32_t gridWidth)
{
	uint8_t rule90Bits = 0b101;
	SetRuleSquareGrid(gridWidth, 1, boost::dynamic_bitset<uint64_t>(3, rule90Bits), 3, 1);
}

void TransitionMatrix::SetRule90Cyclic(uint32_t gridWidth)
{
	uint8_t rule90Bits = 0b101;
	SetRuleToroidGrid(gridWidth, 1, boost::dynamic_bitset<uint64_t>(3, rule90Bits), 3, 1);
}

void TransitionMatrix::SetRuleSquareGrid(uint32_t gridWidth, uint32_t gridHeight, const boost::dynamic_bitset<uint64_t>& rule, int32_t ruleWidth, int32_t ruleHeight)
{
	mGridWidth  = gridWidth;
	mGridHeight = 1;

	uint32_t matrixSize = mGridWidth * mGridHeight;
	mRows.clear();

	for(uint32_t i = 0; i < matrixSize; i++)
	{
		mRows.push_back(boost::dynamic_bitset<uint64_t>(matrixSize));
	}

	for(int32_t rowIndex = 0; rowIndex < gridWidth * gridHeight; rowIndex++)
	{
		int32_t rowCenterX = rowIndex / mGridHeight;
		int32_t rowCenterY = rowIndex % mGridHeight;

		for(int32_t ruleY = 0; ruleY < ruleHeight; ruleY++)
		{
			int32_t ruleYShifted = ruleY - ruleHeight / 2;
			if(rowCenterY + ruleYShifted < 0 || rowCenterY + ruleYShifted >= gridHeight)
			{
				continue;
			}

			for(uint32_t ruleX = 0; ruleX < ruleWidth; ruleX++)
			{
				int32_t ruleXShifted = ruleX - ruleWidth / 2;
				if(rowCenterX + ruleXShifted < 0 || rowCenterX + ruleXShifted >= gridWidth)
				{
					continue;
				}

				uint32_t columnIndex = (uint32_t)(rowCenterY + ruleYShifted) * gridWidth + (uint32_t)(rowCenterX + ruleXShifted);
				mRows[rowIndex][columnIndex] = rule[ruleY * ruleWidth + ruleX];
			} 
		}
	}
}

void TransitionMatrix::SetRuleToroidGrid(uint32_t gridWidth, uint32_t gridHeight, const boost::dynamic_bitset<uint64_t>& rule, int32_t ruleWidth, int32_t ruleHeight)
{
	mGridWidth  = gridWidth;
	mGridHeight = 1;

	uint32_t matrixSize = mGridWidth * mGridHeight;
	mRows.clear();

	for(uint32_t i = 0; i < matrixSize; i++)
	{
		mRows.push_back(boost::dynamic_bitset<uint64_t>(matrixSize));
	}

	for(int32_t rowIndex = 0; rowIndex < gridWidth * gridHeight; rowIndex++)
	{
		int32_t rowCenterX = rowIndex / mGridHeight;
		int32_t rowCenterY = rowIndex % mGridHeight;

		for(int32_t ruleY = 0; ruleY < ruleHeight; ruleY++)
		{
			int32_t ruleYShifted = ruleY - ruleHeight / 2;
			int32_t rowYCycled  = (rowCenterY + ruleYShifted + gridHeight) % gridHeight;

			for(uint32_t ruleX = 0; ruleX < ruleWidth; ruleX++)
			{
				int32_t ruleXShifted = ruleX - ruleWidth / 2;
				int32_t rowXCycled  = (rowCenterX + ruleXShifted + gridWidth) % gridWidth;

				uint32_t columnIndex = (uint32_t)rowYCycled * gridWidth + (uint32_t)rowXCycled;
				if(rule[ruleY * ruleWidth + ruleX])
				{
					mRows[rowIndex][columnIndex].flip();
				}
			} 
		}
	}
}