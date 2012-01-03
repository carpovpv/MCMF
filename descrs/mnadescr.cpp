#include "mnadescr.h"
#include <openbabel/obconversion.h>

MnaDescr::MnaDescr() : DescriptorFactory("MnaDescr")
{

}

const std::vector < double > & MnaDescr::getDescriptors(OBMol * mol, Mode regime)
{
    static int y = 0;
    //const bool regime = true;
    if(regime == Training)
    {
        std::map < OBMol *, std::vector < double > >::iterator it;
        it = m_descrs.find(mol);

        if(it != m_descrs.end())
        {
            return m_descrs[mol];
        }

    }

    std::stringstream st;

    OBConversion conv(NULL, &st);
    conv.AddOption("d", OBConversion::GENOPTIONS);
    conv.SetOutFormat("MNA");

    conv.Write(mol);
    st.seekg(0);

    descrs.clear();

    //названия дескрипторов МНА для текущей структуры.
    std::vector< std::string> curnames;

    while(!st.eof())
    {
        std::string d ;
        std::getline(st, d);
        if(d.length() && d[0] == '-')
            curnames.push_back(d);
    }

    //здесь надо делать ремэппинг....
    //первая структура уже обработана, поэтому необходимо подстройка векторов

    const int n = descrnames.size();
    const int n1 = curnames.size();

    if(m_descrs.size())
    {
        descrs.resize(n, 0.0);

        for(int i =0; i< n1; ++i)
        {
            bool ok = true;
            for(int j =0; j< n; ++j)
            {
                if(curnames[i] == descrnames[j])
                {
                    descrs[j] = 1.0;
                    ok = false;
                    break;
                }
            }
            if(ok)
            {
                //не нашли такой дескриптор, добавляем его
                descrnames.push_back(curnames[i]);
                descrs.push_back(1);
            }
        }

        //увеличиваем вектора предыдущих соединений
        const int newsize = descrnames.size();
        std::map < OBMol * , std::vector< double> >::iterator
        it = m_descrs.begin();

        while(it != m_descrs.end())
        {
            (*it).second.resize(newsize, 0);
            it++;
        }
    }
    else //первый вектор дескрипторов
    {
        descrs.resize(n1, 1.0);
        descrnames = curnames;
        std::cout << descrnames.size() << std::endl;
    }

    if(regime == Training)
        m_descrs[mol] = descrs;

    return descrs;

}
