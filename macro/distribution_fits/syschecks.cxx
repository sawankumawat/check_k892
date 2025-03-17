using namespace std;

void syschecks()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/4rBw_fits/backup2/";
    ifstream file;

    vector<string> names = {"default", "varA1", "varA2", "varA3", "varA4", "varB1", "varB2", "varC1", "varC2", "varC3", "varC4", "varD1", "varD2", "varD3", "varD4", "varTrB1", "varTrB2", "varTrC1", "varTrC2", "varTrD1", "varTrD2", "varToA1", "varToA2", "varToB1", "varToB2", "varToC1", "varToC2", "varToD1", "varToD2", "varToE1", "varToE2", "varToF1", "varToF2"}; // All sources

    // Reading files
    for (int ivars = 0; ivars < names.size(); ivars++)
    {
        file.open(path + Form("fit_params_%s.txt", names[ivars].c_str()));
        if (!file.is_open())
        {
            cerr << "Error opening file: " << names[ivars] << endl;
            return;
        }

        string factor, sign, nameall;
        double value1, value2, uncertainty;
        vector<double> values, uncertainties;
        vector<double> other_values;
        while (getline(file, nameall))
        {
            istringstream iss1(nameall);
            istringstream iss2(nameall);
            if (iss1 >> factor >> value1)
            {
                other_values.push_back(value1);
                // cout<<factor<<" is "<<value1<<endl;
            }
            else if (iss2 >> value1 >> sign >> uncertainty)
            {
                values.push_back(value1);
                uncertainties.push_back(uncertainty);
                // cout<<value1<<" ± "<<uncertainty<<endl;
            }
        }
        string othervalues_names[] = {"significance", "statSignificance", "chi2ndf"};
        cout << names[ivars] << endl;

        if (other_values.size() == 3)
        {
            cout << othervalues_names[0] << " is " << Form("%.2f", other_values[0]) << endl;
            cout << othervalues_names[1] << " is " << Form("%.2f", other_values[1]) << endl;
            cout << othervalues_names[2] << " is " << Form("%.2f", other_values[2]) << endl;
        }
        if (values.size() == 3)
        {
            cout << "fit mass " << Form("%.2f", values[1] * 1000) << " ± " << Form("%.2f", uncertainties[1] * 1000) << endl;
            cout << "fit width " << Form("%.2f", values[2] * 1000) << " ± " << Form("%.2f", uncertainties[2] * 1000) << endl;
        }
        cout << endl;
        file.close();
    }
}