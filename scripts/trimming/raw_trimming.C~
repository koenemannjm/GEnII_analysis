//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 4, 2025   
//
//      
//   The purpose of this script is to trim down and the         
//   raw files for GEnII experiment and produce a
//   single root file that can then be analyzed.
//                                                    
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////     
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TRegexp.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>


void raw_trimming(std::string config, std::string target, std::string pass) {
    // Validate config
    std::vector<std::string> valid_configs = {"kin2", "kin3", "kin4a", "kin4b"};
    if (std::find(valid_configs.begin(), valid_configs.end(), config) == valid_configs.end()) {
        std::cerr << "Error: Invalid configuration '" << config << "'.\n";
        std::cerr << "Allowed configurations are: kin2, kin3, kin4a, kin4b.\n";
        return;
    }

    // Validate target
    std::vector<std::string> valid_targets = {"He3", "H"};
    if (std::find(valid_targets.begin(), valid_targets.end(), target) == valid_targets.end()) {
        std::cerr << "Error: Invalid target '" << target << "'.\n";
        std::cerr << "Allowed configurations are: He3, H.\n";
        return;
    }

    // Validate target
    std::vector<std::string> valid_passes = {"pass1", "pass2"};
    if (std::find(valid_passes.begin(), valid_passes.end(), pass) == valid_passes.end()) {
        std::cerr << "Error: Invalid target '" << pass << "'.\n";
        std::cerr << "Allowed configurations are: pass1, pass2.\n";
        return;
    }

    
    // Set experiment name
    std::string exp_name;
    std::vector<int> selected_numbers;
    if (config == "kin2") {
        exp_name = "GEN2";
	selected_numbers = {2033,2034,2035,2036,2037,2038,2039,2040,2042,2045,2046,2047,2053,2054,2062,2080,2081,2082,2083,2084,2085,2125,2126,2127,2128,2131,2132,2133,2134,2136,2137,2138,2141,2143,2145,2146,2149,2164,2165,2168,2169,2170,2171,2172,2173,2174,2175,2176,2199,2200,2201,2203,2204,2206,2208,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2237,2238,2239,2240,2241,2242,2243,2244,2245,2246,2247,2248,2250,2251,2252,2253,2254,2255,2256,2257,2258,2259,2260,2275,2276,2302,2303,2304,2305,2306,2307,2309,2310,2313,2314,2315,2316,2318,2319,2320,2321,2322,2323};
    }
    else if (config == "kin3") {
         exp_name = "GEN3";
	 selected_numbers = {2506,2507,2511,2512,2516,2517,2518,2519,2520,2521,2522,2523,2524,2525,2526,2527,2528,2529,2531,2535,2536,2537,2539,2540,2541,2542,2543,2544,2545,2551,2552,2553,2554,2555,2556,2557,2558,2559,2560,2561,2562,2563,2564,2565,2566,2568,2569,2574,2570,2571,2572,2574,2575,2576,2577,2578,2579,2580,2581,2583,2584,2585,2586,2587,2613,2614,2615,2616,2617,2618,2619,2620,2621,2622,2623,2624,2625,2626,2627,2628,2629,2630,2631,2632,2633,2634,2635,2637,2638,2646,2647,2648,2649,2650,2651,2658,2659,2660,2661,2662,2667,2668,2669,2670,2671,2674,2676,2677,2678,2679,2680,2681,2683,2684,2685,2687,2688,2689,2690,2691,2692,2693,2694,2695,2696,2697,2698,2699,2700,2701,2702,2703,2704,2705,2706,2707,2709,2710,2711,2712,2713,2714,2715,2716,2717,2718,2719,2720,2721,2722,2723,2724,2725,2726,2727,2728,2729,2730,2763,2764,2765,2766,2767,2768,2769,2770,2959,2964,2966,2967,2968,2970,2971,2973,2974,2976,2981,2983,2984,2985,2986,2988,2989,2990,2995,2996,2997,2998,2999,3000,3001,3003,3004,3005,3006,3007,3008,3009,3011,3013,3014,3016,3017,3018,3019,3020,3021,3022,3024,3025,3026,3027,3028,3029,3030,3035,3036,3037,3038,3039,3043,3045,3046,3047,3048,3049,3050,3051,3052,3054,3055,3058,3059,3061,3062,3063,3064,3065,3068,3069,3075,3076,3077,3079,3080,3082,3089,3090,3091,3092,3093,3108,3109,3110,3111,3112,3113,3114,3115,3116,3117,3118,3119,3120,3121,3122,3123,3124,3126,3127,3134,3135,3136,3138,3139,3140,3141,3149,3150,3151,3152,3158,3159,3160,3161,3162,3163,3164,3165,3166,3167,3170,3171,3173,3174,3176,3177,3178,3179,3180,3181,3182,3183,3184,3185,3186,3188,3190,3191,3192,3193,3194,3195,3196,3197,3198,3199,3200,3201,3202,3203,3204,3205,3212,3213,3214,3215,3217,3218,3219,3220,3221,3222,3223,3224,3225,3226,3227,3228,3229,3230,3231,3239,3240,3241,3242,3243,3244,3245,3246,3247,3248,3249,3250};
    }
    else if (config == "kin4a") {
         exp_name = "GEN4";
         selected_numbers = {3510,3511,3512,3514,3515,3517,3518,3519,3520,3521,3523,3524,3525,3526,3552,3553,3554,3555,3556,3557,3558,3559,3560,3561,3562,3563,3578,3579,3580,3581,3582,3583,3588,3589,3590,3591,3592,3593,3594,3595,3596,3671,3672,3673,3694,3695,3696,3697,3698,3701,3702,3703,3704,3705,3706,3708,3709,3710,3711,3712,3713,3714,3715,3716,3717,3718,3719,3720,3721,3722,3723,3729,3730,3731,3732,3733,3734,3737,3738,3739,3740,3741,3742,3747,3748,3749,3750,3752,3753,3754,3755,3756,3758,3759,3761,3762,3763,3764,3765,3766,3767,3782,3783,3784,3786,3787,3788,3789,3790,3791,3792,3794,3795,3796,3810,3811,3812,3813,3815,3822,3823,3824,3825,3826,3827,3828,3829,3830,3832,3833,3835,3836,3837,3838,3839,3840,3841,3842,3843,3844,3846,3847,3848,3849,3850,3851,3852,3853,3854,3855,3856,3857,3858,3859,3860,3861,3864,3866,3868,3869,3878,3880,3881,3882,3883,3884,3885,3867,3889,3890,3923,3924,3926,3927,3928,3929,3930,3931,3932,3933,3934,3935,3936,3937,3938,3939,3940,3941,3942,3944,3945,3946,3947,3948,3951,3952,3954,3955,3956,3957,3958,3959,3960,3961,3962,3963,3964,3965,3966,3967,3968,3970,3972,3973,3974,3975,3976,3977,3978,3980,3981,3983,3985,3986,3996,3998,4023,4024,4025,4027,4028,4029,4030,4031,4032,4034,4035,4036,4037,4038,4043,4044,4046,4048,4051,4052,4053,4054,4055,4056,4057,4058,4059,4060,4061,4062,4063,4064,4065,4066,4067,4068,4069,4070,4071,4072,4073,4074,4075,4076,4077,4078,4079,4080,4081,4082,4083,4084,4085,4086,4087,4088,4089,4091,4092,4093,4094,4095,4096,4097,4098,4099,4100,4101,4102,4103,4104,4105,4106,4107,4114,4115,4116,4117,4118,4119,4120,4121,4122,4155,4156,4159,4160,4164,4165,4167,4168,4169,4170,4171,4172,4173,4174,4175,4176,4178,4179,4180,4181,4182,4183,4186,4187,4188,4189,4190,4192,4194,4195,4196,4197,4198,4199,4200,4201,4202,4203,4204,4214,4215,4216,4217,4218,4219,4220,4223,4224,4230,4231,4233,4234,4236,4251,4252,4253,4254,4256,4257,4258,4260,4261,4262,4263,4264,4265,4266,4273,4274,4275,4284,4285,4286,4287,4288,4289,4290,4292,4293,4294,4297,4298,4299,4300,4302,4303,4304,4305,4307,4308,4309,4310,4311,4315,4316,4320,4321,4322,4323,4324,4325,4328,4329,4330,4331,4332,4333,4334,4335,4339,4341,4343,4344,4345,4346,4347,4354,4355,4356,4357,4359,4360,4361,4362,4363,4364,4365,4367,4368,4370,4371,4372,4374,4375,4376,4377,4378,4379,4470,4471,4474,4475,4476,4478,4479,4482,4483,4485,4486,4487,4488,4490,4491,4492,4493,4494,4495,4497,4498,4503,4507,4510,4512,4514,4522,4523,4524,4526,4540,4543,4547,4549,4552,4554,4555,4556,4557,4559,4560,4562,4569,4570,4571,4572,4573,4574,4575,4576,4577,4578,4579,4580,4581,4582,4583,4584,4585,4586,4587};
    }
    else if (config == "kin4b") {
         exp_name = "GEN4b";
	 selected_numbers = {5044,5045,5046,5047,5048,5049,5053,5054,5086,5090,5091,5093,5094,5095,5096,5099,5100,5101,5102,5103,5104,5105,5106,5107,5108,5109,5110,5112,5113,5114,5115,5116,5120,5121,5122,5124,5126,5127,5129,5135,5140,5142,5143,5145,5146,5147,5148,5152,5153,5154,5155,5156,5164,5165,5167,5176,5183,5185,5186,5187,5188,5301,5302,5303,5315,5316,5317,5318,5320,5321,5322,5323,5324,5364,5378,5379,5380,5385,5386,5387,5388,5406,5407,5408,5409,5411,5412,5414,5416,5437,5438,5444,5445,5446,5448,5449,5453,5454,5455,5456,5457,5458,5459,5460,5461,5462,5463,5464,5465,5466,5467,5468,5469,5470,5471,5472,5476,5477,5478,5479,5480,5481,5482,5483,5484,5485,5486,5494,5495,5496,5497,5499,5501,5502,5503,5504,5505,5506,5507,5508,5511,5512,5513,5514,5515,5517,5518,5519,5523,5524,5525,5526,5527,5529,5531,5532,5533,5537,5538,5539,5541,5543,5544,5545,5548,5549,5550,5571,5572,5573,5574,5575,5576,5577,5578,5579,5580,5581,5582,5583,5584,5586,5589,5592,5593,5596,5597,5598,5599,5600,5601,5602,5603,5604,5605,5606,5608,5609,5611,5612,5613,5614,5615,5616,5617,5619,5620,5621,5622,5623,5624,5625,5626,5627,5628,5629,5630,5632,5633,5634,5654,5655,5659,5660,5661,5663,5664,5665,5666,5667,5670,5671,5672,5673,5674,5675,5676,5678,5679,5683,5685,5686,5688,5690,5691,5692,5693,5694,5696,5697,5698,5699,5700,5702,5703,5704,5705,5706,5707,5708,5709,5710,5711,5712,5713,5714,5715,5717,5718,5721,5722,5723,5724,5725,5726,5727,5728,5730,5732,5733,5734,5735,5736,5737,5739,5740,5741,5742,5743,5744,5745,5748,5751,5752,5753,5754,5755,5756,5757,5761,5762,5763,5764,5765,5768,5770,5771,5773,5778,5779,5781,5783,5784,5785,5805,5806,5807,5808,5809,5810,5811,5815,5816,5818,5819,5820,5821,5828,5829,5830,5831,5832,5837,5838,5839,5840,5841,5842,5843,5844,5852,5854,5856,5858,5859,5860,5861,5862,5863,5864,5865,5866,5868,5869,5870,5871,5872,5873,5874,5875,5877,5879,5880,5881,5882,5883,5884,5885,5886,5887,5888,5889,5890,5891,5892,5893,5894,5895,5896,5897,5898,5899,5900,5901,5902,5903,5904,5912,5913,5914,5915,5916,5917,5918,5919,5923,5924,5925,5926,5927,5931,5932,5933,5934,5935,5936,5937,5938,5939,5941,5942,5943,5944,5945,5947,5948,5949,5950,5951,5952,5965,5966,5967,5968,5971,5974,5975,5976,5979,5980,5981,5985,5986,5987,5988,5989,5990,5991,5994,5995,5996,5997,5998,5999,6001,6002,6003,6004,6005,6006,6007,6008,6009,6010,6012,6013,6014,6016,6017,6018,6019,6020,6024,6025,6026,6027,6028,6029,6030,6031,6032,6033,6034,6035,6036,6037,6038,6039,6040,6041,6042,6043,6044,6048,6050,6051,6052,6053,6060,6061,6062,6063,6064,6065,6066,6068,6069,6070,6071,6072,6073,6074,6075,6076,6077,6078,6079,6081,6082,6083};
    }

    // Define the branch names you want to extract
    std::vector<std::string> branches_to_extract_scalar = {"g.evtime", "HALLA_p", "bb.grinch_tdc.clus.trackindex", "bb.grinch_tdc.clus.size", "bb.sh.e", "bb.sh.x", "bb.sh.y", "bb.sh.atimeblk", "bb.ps.e", "bb.ps.x", "bb.ps.y", "sbs.hcal.e", "sbs.hcal.x", "sbs.hcal.y", "sbs.hcal.rowblk", "sbs.hcal.colblk", "sbs.hcal.idblk", "sbs.hcal.tdctimeblk", "sbs.hcal.atimeblk", "bb.tr.n", "scalhel.hel", "IGL1I00OD16_16", "Lrb.BPMA.x", "Lrb.BPMA.y", "Lrb.Raster2.rawcur.x", "Lrb.Raster2.rawcur.y", "Lrb.Raster.rawcur.x", "Lrb.Raster.rawcur.y", "sbs.hcal.nclus"}; // Modify as needed

    std::vector<std::string> branches_to_extract_array = {"bb.tr.p", "bb.tr.px", "bb.tr.py", "bb.tr.pz", "bb.tr.vx", "bb.tr.vy", "bb.tr.vz", "bb.tr.tg_x", "bb.tr.tg_y", "bb.tr.tg_th", "bb.tr.tg_ph", "bb.tr.r_x", "bb.tr.r_y", "bb.tr.r_th", "bb.tr.r_ph", "bb.tdctrig.tdcelemID", "bb.tdctrig.tdc",  "bb.hodotdc.clus.bar.tdc.meantime", "bb.gem.track.nhits", "bb.etot_over_p"};
  
   // Directory containing the files
    std::string input_dir = "/cache/halla/sbs/prod/GEnII/" + pass + "/" + exp_name + "/He3/rootfiles/";  // cache directory for Raw data
    TSystemDirectory dir("input", input_dir.c_str());
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "No files found in directory: " << input_dir << std::endl;
        return;
    }

    std::string output_dir = "../../data/raw/" + pass + "/"  + config + "_" + target + "/";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    std::string output_filename = output_dir + "QE_data_" + exp_name + "_sbs100p_nucleon_np" + ".root";

    // Create an output ROOT file
    TFile *outputFile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *outputTree = new TTree("Tout", "Merged data from selected ROOT files");

    std::map<std::string, Double_t> outScalars;
    std::map<std::string, Double_t> outArrays;

    for (const auto& branch : branches_to_extract_scalar) {
        outputTree->Branch(branch.c_str(), &outScalars[branch]);
    }

    for (const auto& branch : branches_to_extract_array) {
        outputTree->Branch(branch.c_str(), &outArrays[branch]);
    }

    int totalFiles = 0;
    for (TIter it(files); TSystemFile *f = (TSystemFile*)it(); ) {
        TString name = f->GetName();
        if (!f->IsDirectory() && name.EndsWith(".root")) {
            for (int num : selected_numbers) {
                if (name.Contains(Form("_%d_", num))) {
                    totalFiles++;
                    break;
                }
            }
        }
    }

    int processed = 0;

    // Loop over all files in the directory and match pattern "_####_"
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            for (int num : selected_numbers) {
                TString pattern = Form("_%d_", num);
                if (fname.Contains(pattern)) {
                    std::string full_path = input_dir + fname.Data();
                    processed++;
		    float percent = (100.0f * processed) / totalFiles;
		    std::cout << "\rProcessing: " << std::fixed << std::setprecision(2) << percent << "% complete" << std::flush;

                    TFile *inputFile = TFile::Open(full_path.c_str(), "READ");
                    if (!inputFile || inputFile->IsZombie()) {
                        std::cerr << "Failed to open file: " << full_path << std::endl;
                        continue;
                    }

                    TTree *inputTree = (TTree*)inputFile->Get("T"); // Replace with actual tree name
                    if (!inputTree) {
                        std::cerr << "Tree not found in file: " << full_path << std::endl;
                        inputFile->Close();
                        continue;
                    }

		    TTreeReader reader(inputTree);
		    std::map<std::string, std::unique_ptr<TTreeReaderValue<Double_t>>> scalarReaders;
		    std::map<std::string, std::unique_ptr<TTreeReaderArray<Double_t>>> arrayReaders;
		    
		    for (const auto& branch : branches_to_extract_scalar) {
			scalarReaders[branch] = std::make_unique<TTreeReaderValue<Double_t>>(reader, branch.c_str());
		    }
		    for (const auto& branch : branches_to_extract_array) {
		        arrayReaders[branch] = std::make_unique<TTreeReaderArray<Double_t>>(reader, branch.c_str());
		    }
		    
                    while (reader.Next()) {
                        for (const auto& [name, readerPtr] : scalarReaders) {
                            outScalars[name] = **readerPtr;
                        }
                        for (const auto& [name, readerPtr] : arrayReaders) {
                            outArrays[name] = readerPtr->GetSize() > 0 ? (*readerPtr)[0] : -9999;
                        }
                        outputTree->Fill();
                    }

                    inputFile->Close();
                    break; // Avoid processing the same file multiple times if multiple numbers match
                }
            }
        }
    }

    std::cout << std::endl;

    
    // Write and close the output file
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();

    std::cout << "Merged ROOT file created: " << output_filename << std::endl;
}
