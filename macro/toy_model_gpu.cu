#include <iostream>
#include <cuda.h>
#include <curand_kernel.h>
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"

// constants
__constant__ double m_mother = 1.713;   // f0(1710) mass
__constant__ double m_daughter1 = 0.493; // Ks
__constant__ double m_daughter2 = 0.493; // Ks

// GPU kernel: generate events
__global__ void decayKernel(int nEvents, float *pT_array, float *recPt_array, unsigned long seed) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nEvents) return;

    // random state
    curandState state;
    curand_init(seed, i, 0, &state);

    // --- generate random kinematics ---
    float pT  = 20.0f * curand_uniform(&state);                    // pT in [0,20]
    float phi = (2.0f * M_PI) * (curand_uniform(&state) - 0.5f);   // phi in [-pi,pi]
    float eta = -0.8f + 1.6f * curand_uniform(&state);             // eta in [-0.8,0.8]

    // mother 4-momentum in lab frame (approximate: using TLorentzVector math on CPU side normally)
    float px = pT * cosf(phi);
    float py = pT * sinf(phi);
    float pz = pT * sinhf(eta);
    float E  = sqrtf(px*px + py*py + pz*pz + m_mother*m_mother);

    // store mother pT
    pT_array[i] = pT;

    // --- simple two-body decay in mother rest frame ---
    float M  = m_mother;
    float m1 = m_daughter1;
    float m2 = m_daughter2;

    // momentum of daughters in CM frame
    float p_star = sqrtf((M*M - (m1+m2)*(m1+m2))*(M*M - (m1-m2)*(m1-m2))) / (2*M);

    // pick random isotropic direction
    float costheta = 2.0f*curand_uniform(&state) - 1.0f;
    float sintheta = sqrtf(1.0f - costheta*costheta);
    float phi_decay = 2.0f*M_PI*curand_uniform(&state);

    float px1 = p_star * sintheta * cosf(phi_decay);
    float py1 = p_star * sintheta * sinf(phi_decay);
    float pz1 = p_star * costheta;
    float E1  = sqrtf(p_star*p_star + m1*m1);

    // reconstruct mother from daughters (back in CM → trivial)
    float rec_px = px1 - px1;  // should cancel, placeholder
    float rec_py = py1 - py1;
    float rec_pz = pz1 - pz1;
    float rec_E  = E1 + sqrtf(p_star*p_star + m2*m2);

    float rec_pT = sqrtf(rec_px*rec_px + rec_py*rec_py);
    recPt_array[i] = rec_pT;
}

int main() {
    int nEvents = 1e8;

    // allocate GPU arrays
    float *d_pT, *d_recPt;
    cudaMalloc(&d_pT, nEvents*sizeof(float));
    cudaMalloc(&d_recPt, nEvents*sizeof(float));

    // launch kernel
    int threads = 256;
    int blocks  = (nEvents + threads - 1) / threads;
    decayKernel<<<blocks, threads>>>(nEvents, d_pT, d_recPt, 1234UL);
    cudaDeviceSynchronize();

    // copy back to CPU
    float *h_pT = new float[nEvents];
    float *h_recPt = new float[nEvents];
    cudaMemcpy(h_pT, d_pT, nEvents*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_recPt, d_recPt, nEvents*sizeof(float), cudaMemcpyDeviceToHost);

    // ROOT histograms
    TH1F *h1 = new TH1F("h1","Mother p_{T};p_{T} (GeV/c);Events",150,0,30);
    TH1F *h2 = new TH1F("h2","Reconstructed p_{T};p_{T} (GeV/c);Events",150,0,30);

    for (int i=0;i<nEvents;i++) {
        h1->Fill(h_pT[i]);
        h2->Fill(h_recPt[i]);
    }

    // draw and save
    TCanvas *c1 = new TCanvas("c1","GPU Hist",800,600);
    h1->SetLineColor(kBlue); h1->Draw();
    c1->SaveAs("mother_pT_gpu.png");

    TCanvas *c2 = new TCanvas("c2","GPU Reco",800,600);
    h2->SetLineColor(kRed); h2->Draw();
    c2->SaveAs("reconstructed_pT_gpu.png");

    // cleanup
    delete[] h_pT;
    delete[] h_recPt;
    cudaFree(d_pT);
    cudaFree(d_recPt);

    return 0;
}
