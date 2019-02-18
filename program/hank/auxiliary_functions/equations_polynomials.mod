// Specifies equations for firstOrderDynamics.mod
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Conditional expectation (#equations = nz * nAssets * nShare)
//----------------------------------------------------------------

@#for iShare in 1 : nShare
	
	@#for iz in 1 : nz 
	
		@#for iAssets in 1 : nAssets
		
			// Compute conditional expectation
			# expectation_@{iz}_@{iAssets}_@{iShare} = exp(0
			@#for iPower in 1 : nAssets
				+ expectationCoefficient_@{iz}_@{iPower}_@{iShare} * expectationPoly_@{iAssets}_@{iPower}
			@#endfor
			);
			
			// Compute savings policy
			# assetsPrime_@{iz}_@{iAssets}_@{iShare} = max(bbBar, 
				(((1-ttau) * w * zGrid_@{iz}) ^ (1 + 1/nnu)) 
				* ((expectation_@{iz}_@{iAssets}_@{iShare} / ppsi) ^ (1/nnu))
                + (1 + r) * assetsGrid_@{iAssets} + shareGrid_@{iShare} * d + T
				- 1 / expectation_@{iz}_@{iAssets}_@{iShare});
				
			// Compute next period's consumption
			@#for izPrime in 1 : nz
			
				# expectationPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} = exp(0
				@#for iPower in 1 : nAssets
					+ expectationCoefficient_@{izPrime}_@{iPower}_@{iShare}(+1) * cos((@{iPower} - 1) * acos(min(max(
						2 * ((assetsPrime_@{iz}_@{iAssets}_@{iShare} - assetsMin) / (assetsMax - assetsMin)) - 1,-1),1)))                                              
				@#endfor
				);
				
				# assetsPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} = max(bbBar, 
				    (((1-ttau) * w(+1) * zGrid_@{izPrime}) ^ (1 + 1/nnu)) 
					* ((expectationPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} / ppsi) ^ (1/nnu))
					+ (1 + r(+1)) * assetsPrime_@{iz}_@{iAssets}_@{iShare} + shareGrid_@{iShare} * d(+1) + T(+1)
					- 1 / expectationPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare});
					
				# auxPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} = 
					-bbBar + (1 + r(+1)) * assetsPrime_@{iz}_@{iAssets}_@{iShare} + shareGrid_@{iShare} * d(+1) + T(+1);
				# consumptionPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} = 
					1 / expectationPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare}
					* (assetsPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} >= bbBar + 1e-8)
					+ (auxPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} + sqrt((auxPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} ^ 2)
					+ 4 * (((1-ttau) * w(+1) * zGrid_@{izPrime})^2) / ppsi)) / 2
					* (assetsPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare} < bbBar + 1e-8);
					
			@#endfor
					
			// Functional equation
			log(expectation_@{iz}_@{iAssets}_@{iShare}) = log(bbeta * (1 + r(+1)) * (0
			@#for izPrime in 1 : nz
				+ zTransition_@{iz}_@{izPrime} / consumptionPrime_@{izPrime}_@{iz}_@{iAssets}_@{iShare}
			@#endfor
			));
			
		@#endfor

	@#endfor
	
@#endfor

//----------------------------------------------------------------
// Compute various objects over quadrature grid for integrating distribution
//----------------------------------------------------------------

@#for iShare in 1 : nShare

	@#for iz in 1 : nz

		@#for iAssets in 1 : nAssetsQuadrature
		
			// Compute conditional expectation
			# expectationQuadrature_@{iz}_@{iAssets}_@{iShare} = exp(0
			@#for iPower in 1 : nAssets
				+ expectationCoefficient_@{iz}_@{iPower}_@{iShare} * quadraturePoly_@{iAssets}_@{iPower}
			@#endfor
			);
			
			// Compute savings policy
			# assetsPrimeQuadrature_@{iz}_@{iAssets}_@{iShare} = max(bbBar, 
				(((1-ttau) * w * zGrid_@{iz}) ^ (1 + 1/nnu)) 
				* ((expectationQuadrature_@{iz}_@{iAssets}_@{iShare} / ppsi) ^ (1/nnu))
                + (1 + r) * quadratureGrid_@{iAssets} + shareGrid_@{iShare} * d + T
				- 1 / expectationQuadrature_@{iz}_@{iAssets}_@{iShare});
				
			// Compute labor supply (labor*z)
			# aux_@{iz}_@{iAssets}_@{iShare} = -bbBar + (1 + r) * assetsPrimeQuadrature_@{iz}_@{iAssets}_@{iShare} + shareGrid_@{iShare} * d + T;
			# laborQuadrature_@{iz}_@{iAssets}_@{iShare} = (1-ttau) * w * (zGrid_@{iz} ^ 2)
				* expectationQuadrature_@{iz}_@{iAssets}_@{iShare} / ppsi 
				* (assetsPrimeQuadrature_@{iz}_@{iAssets}_@{iShare} >= bbBar + 1e-8)
				+ (-aux_@{iz}_@{iAssets}_@{iShare} + sqrt((aux_@{iz}_@{iAssets}_@{iShare} ^ 2)
				+ 4 * (((1-ttau) * w * zGrid_@{iz})^2) / ppsi))/ (2 * (1-ttau) * w)
				* (assetsPrimeQuadrature_@{iz}_@{iAssets}_@{iShare} < bbBar + 1e-8);
				
			// PDF of distribution
			# measurePDF_@{iz}_@{iAssets}_@{iShare} = exp(0 + measureCoefficient_@{iz}_1_@{iShare} * (quadratureGrid_@{iAssets} - 
				moment_@{iz}_1_@{iShare}(-1))
				@#for iMoment in 2 : nMeasure
					+ measureCoefficient_@{iz}_@{iMoment}_@{iShare} * ((quadratureGrid_@{iAssets} - moment_@{iz}_1_@{iShare}(-1)) ^ @{iMoment} - 
					moment_@{iz}_@{iMoment}_@{iShare}(-1))
				@#endfor
				);
				
		@#endfor
		
		// Total mass of distribution
		# totalMass_@{iz}_@{iShare} = 0
		@#for iAssets in 1 : nAssetsQuadrature
			+ quadratureWeights_@{iAssets} * measurePDF_@{iz}_@{iAssets}_@{iShare}
		@#endfor
		;
		
	@#endfor
	
@#endfor

//----------------------------------------------------------------
// Compute various objects at borrowing constraint for integrating distribution
//----------------------------------------------------------------

@#for iShare in 1 : nShare

	@#for iz in 1 : nz

		// Compute conditional expectation
		# expectationBC_@{iz}_@{iShare} = exp(0
		@#for iPower in 1 : nAssets
			+ expectationCoefficient_@{iz}_@{iPower}_@{iShare} * bcPoly_@{iPower}
		@#endfor
		);
		
		// Compute savings policy
		# assetsPrimeBC_@{iz}_@{iShare} = max(bbBar, 
			(((1-ttau) * w * zGrid_@{iz}) ^ (1 + 1/nnu)) 
			* ((expectationBC_@{iz}_@{iShare} / ppsi) ^ (1/nnu))
			+ (1 + r) * bbBar + shareGrid_@{iShare} * d + T
			- 1 / expectationBC_@{iz}_@{iShare});
			
		// Compute labor supply (labor*z)
		# auxBC_@{iz}_@{iShare} = -bbBar + (1+r) * assetsPrimeBC_@{iz}_@{iShare} + shareGrid_@{iShare} * d + T;
		# laborBC_@{iz}_@{iShare} = (1-ttau) * w * (zGrid_@{iz} ^ 2)
			* expectationBC_@{iz}_@{iShare} / ppsi 
			* (assetsPrimeBC_@{iz}_@{iShare} >= bbBar + 1e-8)
			+ (-auxBC_@{iz}_@{iShare} + sqrt((auxBC_@{iz}_@{iShare} ^ 2) 
			+ 4 * (((1-ttau) * w * zGrid_@{iz})^2) / ppsi))/ (2 * (1-ttau) * w)
			* (assetsPrimeBC_@{iz}_@{iShare} < bbBar + 1e-8);
			
	@#endfor
	
@#endfor

//----------------------------------------------------------------
// Relationship between moments of distribution and parameters
// (#equations = nz * nMeasure * nShare)
//----------------------------------------------------------------

@#for iShare in 1 : nShare

	@#for iz in 1 : nz
		
		// First moments (uncentered)
		moment_@{iz}_1_@{iShare}(-1) = (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * quadratureGrid_@{iAssets} * measurePDF_@{iz}_@{iAssets}_@{iShare}
			@#endfor
			) / totalMass_@{iz}_@{iShare};
		
		// Higher order moments (centered)
		@#for iMoment in 2 : nMeasure
		moment_@{iz}_@{iMoment}_@{iShare}(-1) = (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{iz}_@{iAssets}_@{iShare} * ((quadratureGrid_@{iAssets} - 
					moment_@{iz}_1_@{iShare}(-1)) ^ @{iMoment})
			@#endfor
			) / totalMass_@{iz}_@{iShare};
			
		@#endfor
			
	@#endfor

@#endfor

//----------------------------------------------------------------
// Law of motion for density away from borrowing constraint 
// (#equations = nz * nMeasure * nShare)
//----------------------------------------------------------------

@#for iShare in 1 : nShare

	@#for iz in 1 : nz
		
		// First moment (uncentered)
		moment_@{iz}_1_@{iShare} = (0
		@#for izTilde in 1 : nz
			+ ((1 - mHat_@{izTilde}_@{iShare}(-1)) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz} * (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{izTilde}_@{iAssets}_@{iShare} *
					assetsPrimeQuadrature_@{izTilde}_@{iAssets}_@{iShare}
			@#endfor
			) / totalMass_@{izTilde}_@{iShare}) + mHat_@{izTilde}_@{iShare}(-1) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz} 
				* assetsPrimeBC_@{izTilde}_@{iShare}
		@#endfor
		) / zMass_@{iz};
		
		// Higher order moments (centered)
		@#for iMoment in 2 : nMeasure
			moment_@{iz}_@{iMoment}_@{iShare} = (0
			@#for izTilde in 1 : nz
				+ ((1 - mHat_@{izTilde}_@{iShare}(-1)) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz} * (0
				@#for iAssets in 1 : nAssetsQuadrature
					+ quadratureWeights_@{iAssets} * measurePDF_@{izTilde}_@{iAssets}_@{iShare} * 
						(assetsPrimeQuadrature_@{izTilde}_@{iAssets}_@{iShare} - moment_@{iz}_1_@{iShare}) ^ @{iMoment}
				@#endfor
				) / totalMass_@{izTilde}_@{iShare}) + mHat_@{izTilde}_@{iShare}(-1) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz}
					* (assetsPrimeBC_@{izTilde}_@{iShare} - moment_@{iz}_1_@{iShare}) ^ @{iMoment}
			@#endfor
			) / zMass_@{iz};
		@#endfor

	@#endfor

@#endfor

//----------------------------------------------------------------
// Law of motion for mass at borrowing constraint
// (#equations = nz * nShare)
//----------------------------------------------------------------

@#for iShare in 1 : nShare

	@#for iz in 1 : nz

		mHat_@{iz}_@{iShare} = (0
		@#for izTilde in 1 : nz
			+ ((1 - mHat_@{izTilde}_@{iShare}(-1)) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz} * (0
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{izTilde}_@{iAssets}_@{iShare} *
					(assetsPrimeQuadrature_@{izTilde}_@{iAssets}_@{iShare} <= bbBar + 1e-8)
			@#endfor
			) / totalMass_@{izTilde}_@{iShare}) + mHat_@{izTilde}_@{iShare}(-1) * zMass_@{izTilde} * zTransition_@{izTilde}_@{iz} 
				* (assetsPrimeBC_@{izTilde}_@{iShare} <= bbBar + 1e-8)
		@#endfor
		) / zMass_@{iz};

	@#endfor

@#endfor

//----------------------------------------------------------------
// Firm conditions (# equations = 2)
//----------------------------------------------------------------

d = (1 - w / A - ttheta /2 * (ppi ^ 2)) * A * N;

ppi = (eepsilon - 1) / ttheta * (log(w / w_SS) - log(A / A_SS)) + 1 / (1 + r_SS) * ppi(+1);

//----------------------------------------------------------------
// Market clearing conditions (# equations = 2)
//----------------------------------------------------------------

N = 0
	@#for iShare in 1 : nShare
		+ shareMass_@{iShare}*(
		@#for iz in 1 : nz
			+ zMass_@{iz} * ((1 - mHat_@{iz}_@{iShare}(-1))*(
			@#for iAssets in 1 : nAssetsQuadrature
				+ quadratureWeights_@{iAssets} * measurePDF_@{iz}_@{iAssets}_@{iShare} * laborQuadrature_@{iz}_@{iAssets}_@{iShare}
			@#endfor
			) / totalMass_@{iz}_@{iShare} + mHat_@{iz}_@{iShare}(-1) * laborBC_@{iz}_@{iShare}) 
		@#endfor
		)
	@#endfor
	;

B_SS + 
	@#for iShare in 1 : nShare
		+ shareMass_@{iShare}*(
		@#for iz in 1 : nz
			+ zMass_@{iz} * ((1 - mHat_@{iz}_@{iShare}(-1)) * moment_@{iz}_1_@{iShare}(-1)
			+ mHat_@{iz}_@{iShare}(-1) * bbBar) 
		@#endfor
		)
	@#endfor
	= 0;

//----------------------------------------------------------------
// Policy equations (# equations = 2)
//----------------------------------------------------------------

// Fiscal policy
T = r * B_SS - G_SS + ttau * w * N;

// Monetary policy
i = r_SS + pphi * ppi + epsilonM;

// Fisher equation
1 + r = (1 + i(-1)) / (1 + ppi);

//----------------------------------------------------------------
// Driving process equations (# equations = 1)
//----------------------------------------------------------------

// Law of motion for aggregate TFP
log(A) = rrhoTFP * log(A(-1)) + (1-rrhoTFP) * log(A_SS) + ssigmaTFP * epsilonA;

