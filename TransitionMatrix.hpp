#pragma once

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using StateVector = boost::dynamic_bitset<uint64_t>;

class TransitionMatrix
{
public:
	TransitionMatrix();
	~TransitionMatrix();

	void SetIdentity(uint32_t gridWidth, uint32_t gridHeight);

	void SetRule150Square(uint32_t gridWidth);
	void SetRule150Cyclic(uint32_t gridWidth);

	TransitionMatrix Mul(const TransitionMatrix& right);
	StateVector MulState(const StateVector& state);

	TransitionMatrix CalcMatrixPower(const boost::multiprecision::cpp_int& matrixPower);

private:
	void SetRuleSquareGrid(uint32_t gridWidth, uint32_t gridHeight, const boost::dynamic_bitset<uint64_t>& rule, int32_t ruleWidth, int32_t ruleHeight);
	void SetRuleToroidGrid(uint32_t gridWidth, uint32_t gridHeight, const boost::dynamic_bitset<uint64_t>& rule, int32_t ruleWidth, int32_t ruleHeight);

private:
	uint32_t mGridWidth  = 0;
	uint32_t mGridHeight = 0;

	std::vector<boost::dynamic_bitset<uint64_t>> mRows;
};