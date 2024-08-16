#include <iostream>
#include "TransitionMatrix.hpp"
#include <bit>
#include <regex>
#include <optional>
#include <format>

template<> struct std::formatter<boost::multiprecision::cpp_int>: std::formatter<std::string> 
{
	auto format(const boost::multiprecision::cpp_int& val, std::format_context& context) const
	{
		return std::formatter<std::string>::format(val.str(), context);
	}
};

enum class LaunchMode
{
	CalcCycleOffset,
	VerifyCycleLength,
};

enum class GridTopology
{
	Square,
	Torus,
};

enum class RuleType
{
	Rule150
};

struct LaunchOptions
{
	LaunchMode LaunchMode = LaunchMode::CalcCycleOffset;

	boost::multiprecision::cpp_int CycleLengthToVerify = 0;
	boost::multiprecision::cpp_int CycleLengthOffset   = 0;

	RuleType Rule = RuleType::Rule150;

	uint32_t GridWidth  = 5;
	uint32_t GridHeight = 5;

	GridTopology Topology = GridTopology::Square;
};

std::optional<LaunchOptions> ParseCommandLineArgs(int argc, char* argv[])
{
	LaunchOptions result;

	bool calcCycleOffset  = false;

	int currArg = 1;
	while(currArg < argc)
	{
		std::string_view currArgStr = argv[currArg];

		if(currArgStr == "--size")
		{
			currArg++;
			
			if(currArg >= argc)
			{
				std::cout << "Please enter grid size!" << std::endl;
				return std::nullopt;
			}
			else
			{
				std::string currStr(argv[currArg]);

				std::regex singleDigitSizeRegex("(\\d+)");
				std::regex doubleDigitSizeRegex("(\\d+)x(\\d+)");
				
				std::smatch singleDigitMatch;
				std::smatch doubleDigitMatch;

				try
				{
					if(std::regex_match(currStr, singleDigitMatch, singleDigitSizeRegex))
					{
						//The entered size is a single number. Parse as a 1D grid
						uint32_t width = std::stoul(singleDigitMatch[1].str());

						result.GridWidth  = width;
						result.GridHeight = 1;
					}
					else if(std::regex_match(currStr, doubleDigitMatch, doubleDigitSizeRegex))
					{
						//The entered size is two numbers. Parse as an n x m grid
						result.GridWidth  = std::stoul(doubleDigitMatch[1].str());
						result.GridHeight = std::stoul(doubleDigitMatch[2].str());
					}
				}
				catch(...)
				{
					std::cout << "Please enter valid board size!" << std::endl;
					return std::nullopt;
				}
			}
		}

		else if(currArgStr == "--offset")
		{
			currArg++;
			
			if(currArg >= argc)
			{
				std::cout << "Please enter cycle length offset!" << std::endl;
				return std::nullopt;
			}
			else
			{
				try
				{
					result.CycleLengthOffset = boost::multiprecision::cpp_int(argv[currArg]);
				}
				catch(std::runtime_error e)
				{
					std::cout << "Please enter valid cycle length offset!" << std::endl;
					return std::nullopt;
				}
			}
		}

		else if(currArgStr == "--rule-150")
		{
			result.Rule = RuleType::Rule150;
		}

		else if(currArgStr == "--topology")
		{
			currArg++;

			if(currArg >= argc)
			{
				std::cout << "Please enter topology!" << std::endl;
				return std::nullopt;
			}
			else if(strcmp(argv[currArg], "square") == 0)
			{
				result.Topology = GridTopology::Square;
			}
			else if(strcmp(argv[currArg], "torus") == 0)
			{
				result.Topology = GridTopology::Torus;
			}
			else
			{
				std::cout << "Unknown topology: " << argv[currArg] << std::endl;
				return std::nullopt;
			}
		}

		else if(currArgStr == "--calc_cycle_offset")
		{
			currArg++;

			calcCycleOffset = true;
			
			if(currArg >= argc)
			{
				std::cout << "Please enter a conservative cycle length!" << std::endl;
				return std::nullopt;
			}
			else
			{
				try
				{
					result.CycleLengthToVerify = boost::multiprecision::cpp_int(argv[currArg]);
				}
				catch(std::runtime_error e)
				{
					std::cout << "Please enter a valid conservative cycle length!" << std::endl;
					return std::nullopt;
				}
			}
		}

		else if(currArgStr == "--verify_cycle_length")
		{
			currArg++;
			
			if(currArg >= argc)
			{
				std::cout << "Please enter cycle length to verify!" << std::endl;
				return std::nullopt;
			}
			else
			{
				try
				{
					result.CycleLengthToVerify = boost::multiprecision::cpp_int(argv[currArg]);
				}
				catch(std::runtime_error e)
				{
					std::cout << "Please enter valid cycle length!" << std::endl;
					return std::nullopt;
				}

				if(result.CycleLengthToVerify <= 0)
				{
					std::cout << "The cycle length must be a positive integer!" << std::endl;
					return std::nullopt;	
				}
			}
		}

		else
		{
			std::cout << "Unknown option: " << argv[currArg] << std::endl;
			return std::nullopt;
		}

		currArg++;
	}

	if(calcCycleOffset)
	{
		result.LaunchMode = LaunchMode::CalcCycleOffset;
	}
	else if(result.CycleLengthToVerify != 0)
	{
		result.LaunchMode = LaunchMode::VerifyCycleLength;
	}

	return result;
}

bool VerifyCycleLength(uint32_t gridWidth, uint32_t gridHeight, GridTopology topology, RuleType rule, const boost::multiprecision::cpp_int& cycleLength, const boost::multiprecision::cpp_int& cycleLengthOffset)
{
	TransitionMatrix mat;
	
	if(topology == GridTopology::Square)
	{
		if(rule == RuleType::Rule150)
		{
			mat.SetRule150Square(gridWidth);
		}
	}
	else if(topology == GridTopology::Torus)
	{
		if(rule == RuleType::Rule150)
		{
			mat.SetRule150Cyclic(gridWidth);
		}
	}

	boost::dynamic_bitset<uint64_t> vectorTest(gridWidth * gridHeight, 0);
	vectorTest.set(0, true);

	//Test that the (A^p)*((A^offset)*b) == (A^offset)*b, i.e. this is indeed the solution period
	{
		TransitionMatrix offsetMatrix = mat.CalcMatrixPower(cycleLengthOffset);
		boost::dynamic_bitset<uint64_t> vectorResBase = offsetMatrix.MulState(vectorTest);

		TransitionMatrix cycleLengthMatrix = mat.CalcMatrixPower(cycleLength + cycleLengthOffset);
		boost::dynamic_bitset<uint64_t> vectorRes = cycleLengthMatrix.MulState(vectorTest);

		if(vectorRes != vectorResBase)
		{
			return false;
		}
	}

	return true;
}

int main(int argc, char *argv[])
{
	auto launchOptionsOpt = ParseCommandLineArgs(argc, argv);
	if(!launchOptionsOpt.has_value())
	{
		return 1;
	}

	auto launchOptions = launchOptionsOpt.value();

	if(launchOptions.LaunchMode == LaunchMode::VerifyCycleLength)
	{
		uint32_t gridWidth  = launchOptions.GridWidth;
		uint32_t gridHeight = launchOptions.GridHeight;

		if(!VerifyCycleLength(gridWidth, gridHeight, launchOptions.Topology, launchOptions.Rule, launchOptions.CycleLengthToVerify, launchOptions.CycleLengthOffset))
		{
			std::cout << "Cycle length verification error!" << std::endl;
		}
		else
		{
			std::cout << "Cycle length verification success!" << std::endl;
		}
	}
	else if(launchOptions.LaunchMode == LaunchMode::CalcCycleOffset)
	{
		uint32_t gridWidth = launchOptions.GridWidth;
		uint32_t gridHeight = launchOptions.GridHeight; //Only square boards are supported by this

		if(VerifyCycleLength(gridWidth, gridHeight, launchOptions.Topology, launchOptions.Rule, launchOptions.CycleLengthToVerify, 0))
		{
			std::cout << std::format("Cycle length offset is 0.", gridWidth, gridHeight) << std::endl;
		}
		else
		{
			boost::multiprecision::cpp_int maxOffsetToTry = 1;
			while(true)
			{
				if(VerifyCycleLength(gridWidth, gridHeight, launchOptions.Topology, launchOptions.Rule, launchOptions.CycleLengthToVerify, maxOffsetToTry))
				{
					break;
				}

				maxOffsetToTry *= 2;
			}

			if(maxOffsetToTry == 1)
			{
				std::cout << std::format("Cycle length offset is 1.", gridWidth, gridHeight) << std::endl;
			}
			else
			{
				boost::multiprecision::cpp_int minOffsetToTry = maxOffsetToTry / 2;
				while(maxOffsetToTry != minOffsetToTry)
				{
					boost::multiprecision::cpp_int currOffsetToTry = (maxOffsetToTry + minOffsetToTry) / 2;
					if(VerifyCycleLength(gridWidth, gridHeight, launchOptions.Topology, launchOptions.Rule, launchOptions.CycleLengthToVerify, currOffsetToTry))
					{
						maxOffsetToTry = currOffsetToTry;
					}
					else
					{
						minOffsetToTry = currOffsetToTry + 1;
					}
				}

				std::cout << std::format("Cycle length offset is {}.", maxOffsetToTry) << std::endl;
			}
		}
	}

	return 0;
}
