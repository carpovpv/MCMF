#ifndef PARSER_H
#define PARSER_H

//Структура с нашими параметрами.

#define D_ELECTROSTATIC 1
#define D_STERIC        2
#define D_HYDROPHOBIC   3
#define D_GAUSS         4
#define D_FP2           5
#define D_MNA           6
#define D_SPECTROPHORES 7
#define D_LINEAR        8
#define D_TANIMOTO      9
#define D_UNKNOWNDESCR  10
#define D_UNKNOWNKERNEL 11
#define D_ABRAHAMA      12
#define D_ABRAHAMB      13
#define D_ABRAHAME      14
#define D_ABRAHAMS      15
#define D_HYDROPHOBICV  16
#define D_STERICK       17

struct ConditionKernel
{
        int kernel;    //ядро
        int descr;     //его дескрипторы
};

struct ConditionParams
{
        std::vector< double> params;				//параметры для линейной комбинации
        std::vector< struct ConditionKernel > kernels;		//вектор ядер
        std::string sdf_test;					//имя файла с тестовыми структурами
        std::string sdf_train;					//имя файла с учебными структурами
        std::string results;					//имя файла для записи результатов
        std::string model;					//имя файла для сохранения модели
        int cv;							//количество выборок для перекрестного контроля
        int max_iter;						//максимальное количество рандомных вариантов обучения
        std::string machine;					//имя машины
        bool prognosis;						//режим прогнозатора
        bool help;						//помощь
};

extern "C" bool parse_command_line(const char * str);
extern struct ConditionParams cond;

#endif // PARSER_H
