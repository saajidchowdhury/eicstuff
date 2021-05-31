R__LOAD_LIBRARY(libeicsmear);
void read() {
    double alpha = 1.0/137.036, 
           scm = 4.0*10*100,
           fbgev = 1.0/(0.3894e12),
           mp = 0.9383;

    double q2min = 1e-1, q2max = 1e4;
    int nq2 = 25; 
    double logbwq2 = (log10(q2max) - log10(q2min)) / nbinsq2;
    double q2low[nq2], q2hi[nq2];
    for (int i = 0; i < nq2; i++) {
        double loglowq2 = log10(q2min) + i*logbwq2;
        double loghiq2 = log10(q2min) + (i+1)*logbwq2;
        q2low[i] = pow(10, lowlowq2);
        q2hi[i] = pow(10, lowhiq2);
    }

    double xmin = 1e-5, xmax = 1;
    int n = 25;
    double logbw = (log10(xmax) - log10(xmin)) / nbins;
    double xlow[n], xhi[n];
    for (int i = 0; i < n; i++) {
        double loglow = log10(xmin) + i*logbw;
        double loghi = log10(xmin) + (i+1)*logbw;
        xlow[i] = pow(10, loglow);
        xhi[i] = pow(10, loghi);
    }

    double yield1[n][nq2], rcs1[n][nq2], lum1[n][nq2], apv1[n][nq2],
           yield2[n][nq2], rcs2[n][nq2], lum2[n][nq2], apv2[n][nq2],
           y[n][nq2],
           q2mid[nq2], xmid[n],
           q2w[nq2], xw[n],
           rcsf[n][nq2];
    bool cutul[n][nq2],
         cutlr[n][nq2];
    double ymax = 1, ymin = 0.001, w2min = 10;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nq2; j++) {
            yield1[i][j] = 0; rcs1[i][j] = 0; lum1[i][j] = 0; apv1[i][j] = 0;
            yield2[i][j] = 0; rcs2[i][j] = 0; lum2[i][j] = 0; apv2[i][j] = 0;

            q2mid[j] = (q2hi[j] + q2low[j]) / 2;
            xmid[i] = (xhi[i] + xlow[i]) / 2;
            q2w[j] = q2hi[j] - q2low[j];
            xw[i] = xhi[i] - xlow[i];

            y[i][j] = q2mid[j] / (xmid[i]*scm);
            rcsf[i][j] = fbgev * pow(q2mid[j], 2) * xmid[i] / 
                         (TMath::TwoPi() * pow(alpha, 2) * (1 + pow(1-y[i][j], 2)));

            double ytemp = q2hi[j] / (xlow[i]*scm);
            if (ytemp > ymax) {
                cutul[i][j] = 0;
            } else {
                cutul[i][j] = 1;
            }
            ytemp = q2low / (xhi[i]*scm);
            double w2temp = mp*mp + q2low[j] * (1 / xhi[i] - 1);
            if (ytemp < ymin or w2temp < w2min) {
                cutlr[i][j] = 0;
            } else {
                cutlr[i][j] = 1;
            }
        }
    }

    int nplot = 8;
    int q2binplot[nplot] = {6, 8, 9, 10, 11, 13, 15, 18};

    erhic::EventDjangoh *event1(NULL);
    erhic::ParticleMC *particle1(NULL);
    TChain *t1 = new TChain("EICTree");
    t1->Add("/gpfs/mnt/gpfs02/eic/saajidchowdhury/outfiles/pol-1/djangoh.NC.Apve.noRad.20x250_evt.root");
    t1->SetBranchAddress("event", &event1);
    t1->GetEntry(0);
    double Q2 = event1->GetQ2();
    double x = event1->GetX();
    double y = event1->GetY();
    double W2 = event1->GetW2();
    cout << Q2 << " " << x << " " << y << " " << W2 << endl;

    erhic::EventDjangoh *event2(NULL);
    erhic::ParticleMC *particle2(NULL);
    TChain *t2 = new TChain("EICTree");
    t2->Add("/gpfs/mnt/gpfs02/eic/saajidchowdhury/outfiles/pol+1/djangoh.NC.Apve.noRad.20x250_evt.root");
    t2->SetBranchAddress("event", &event2);
    t2->GetEntry(0);
    Q2 = event2->GetQ2();
    x = event2->GetX();
    y = event2->GetY();
    W2 = event2->GetW2;
    cout << Q2 << " " << x << " " << y << " " << W2 << endl;
}
