#include "mnadescr.h"
#include <openbabel/obconversion.h>

MnaDescr::MnaDescr() : DescriptorFactory("MnaDescr")
{
    descrcode = D_MNA;
}

void MnaDescr::load(FILE *fp)
{
    int ndescr;
    fscanf(fp, "%d", &ndescr);
    char name[255];

    descrnames.resize(ndescr);
    descrs.resize(ndescr, 0);

    for(int i=0; i< ndescr; ++i)
    {
        fscanf(fp, "%s", name);
        descrnames[i] = name;
        //printf("Loading: %s\n", name);
    }
    fgets(name, 255, fp);
    fprintf(stderr,"Loaded %d descriptors.\n", descrnames.size());
}

const std::vector < double > & MnaDescr::getDescriptors(OBMol * mol, Mode regime)
{
    if(regime == Training) //если обучение, то часть векторов кешируется
    {
        std::map < OBMol *, std::vector < double > >::iterator it;
        it = m_descrs.find(mol);

        if(it != m_descrs.end())
            return m_descrs[mol];
    }
    else
    {
        OBGenericData * data = mol->GetData("prognosis");
        if(data)
        {
            std::string c = data->GetValue();
            long cur = atol(c.c_str());

            if( prev == cur)
                return descrs;
            prev = cur;
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
        {
            if(find(curnames.begin(), curnames.end(), d) == curnames.end())
                curnames.push_back(d);
        }
    }

    //здесь надо делать ремэппинг....
    //первая структура уже обработана, поэтому необходимо подстройка векторов
    //если мы в режиме прогноза, то не нужно расширение вектора, потому что он
    //заполнен при обучении

    const int n = descrnames.size();
    const int n1 = curnames.size();

    if(regime == Training)
    {
        //if(!prognosis)
        {
            if(m_descrs.size() || prognosis)
            {
                descrs.assign(n, 0.0);

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
                if(!prognosis)
                {
                    descrs.assign(n1, 1.0);
                    descrnames = curnames;
                }
            }
            m_descrs[mol] = descrs;
        }
    }
    else //прогнозатор
    {
        descrs.assign(n, 0.0);
        for(int i =0; i< n; ++i)
        {
            for(int j =0; j< n1; ++j)
            {
                if(curnames[j] == descrnames[i])
                {
                    descrs[i] = 1.0;
                    break;
                }
            }
        }

        //for(int i=0; i< n; ++i)
        //    std::cout << descrs[i] << " ";
        //std::cout << std::endl;
    }

    return descrs;
}
