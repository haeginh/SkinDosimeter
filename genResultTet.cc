#include "fstream"
#include "iostream"
#include "vector"
#include "algorithm"
#include "map"
#include "sstream"
#include "cmath"

using namespace std;
int main(int argc,char** argv)
{
    if(argc!=3)
    {
        cout<<"USAGE: ./genResultTet [tet name] [result file]"<<endl;
        return 1;
    }

    string tetName = argv[1];
    string resultFile = argv[2];
    int num;
    string tmp;

    cout<<"reading ele file..."<<flush;
    vector<vector<int>> ELE;
    vector<int> extNode;
    ifstream ifsELE(tetName+".ele");
    ifsELE>>num>>tmp>>tmp;
    for(int i=0;i<num;i++)
    {
        int a, b, c, d, id;
        ifsELE>>tmp>>a>>b>>c>>d>>id;
        if(id>0) continue;
        ELE.push_back({a, b, c, d, (int)floor(-id*0.1)});
        extNode.push_back(a);
        extNode.push_back(b);
        extNode.push_back(c);
        extNode.push_back(d);
    }
    ifsELE.close();
    cout<<"done"<<endl;
    sort(extNode.begin(), extNode.end());
    extNode.erase(unique(extNode.begin(), extNode.end()), extNode.end());
    ifstream ifsNODE(tetName+".node");
    ifsNODE>>num>>tmp>>tmp>>tmp;
    cout<<"# of nodes were decreased: "<<num<<" -> "<<extNode.size()<<endl;
    cout<<"reading node file..."<<flush;
    vector<vector<double>> NODE;
    map<int, int> whole2ext;
    for(int i=0, n=0;i<num;i++)
    {
        double x, y, z;
        ifsNODE>>tmp>>x>>y>>z;
        if(extNode[n]==i)
        {
            NODE.push_back({x,y,z});
            whole2ext[i] = n;
            n++;
        }
    }
    ifsNODE.close();
    cout<<"done"<<endl;

    cout<<"reading result file "<<flush;
    ifstream ifsRslt(resultFile);
    getline(ifsRslt, tmp);
    stringstream ss(tmp);
    ss>>tmp;
    int numSkinDet(-1);
    while(ss>>num>>tmp)
        numSkinDet = numSkinDet<(num+1)? (num+1):numSkinDet;
    cout<<"(skin det: "<<numSkinDet<<")"<<endl;
    map<string, vector<double>> results;
    while(!ifsRslt.eof())
    {
        string name;
        if(!(ifsRslt>>name)) break;
        cout<<"\t"<<name<<": "<<flush;
        double min(__DBL_MAX__), max(0);
        for(int i=0;i<numSkinDet;i++)
        {
            double val, err;
            ifsRslt>>val>>err;
            results[name].push_back(val);
            min = min>val? val:min;
            max = max<val? val:max;
        }
        for(int i=0;i<numSkinDet;i++)
        {
            results[name][i] = (results[name][i]-min)/(max-min)*10000.;
        }
        cout<<min<<" ~ "<<max<<endl;
    }
    ifsRslt.close();
    cout<<"done ("<<results.size()<<")"<<endl;

    cout<<"writing NODE file...("+tetName+"_skin.node)..."<<flush;
    ofstream ofsNODE(tetName+"_skin.node");
    ofsNODE<<NODE.size()<<"   3   0   0"<<endl;
    for(int i=0;i<NODE.size();i++)
    {
        ofsNODE<<i<<" "<<NODE[i][0]<<" "<<NODE[i][1]<<" "<<NODE[i][2]<<endl;
    }
    ofsNODE.close();
    cout<<"done"<<endl;
    for(auto iter:results)
    {
        cout<<"writing "+tetName+"_"+iter.first+".ele..."<<flush;
        ofstream ofsELE(tetName+"_"+iter.first+".ele");
        ofsELE<<ELE.size()<<"   4   1"<<endl;
        for(int i=0;i<ELE.size();i++)
        {
            ofsELE<<i<<" "<<whole2ext[ELE[i][0]]<<" "<<whole2ext[ELE[i][1]]<<" "<<whole2ext[ELE[i][2]]<<" "<<whole2ext[ELE[i][3]]<<" "<<iter.second[ELE[i][4]]<<endl;
        }
        ofsELE.close();
        cout<<"done"<<endl;
    }
    
    return 0;
}