#ifndef PARSER_H
#define PARSER_H
//Структура с нашими параметрами.

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

extern "C" bool parse_command_line(char * str);
extern struct ConditionParams cond;

#endif // PARSER_H
