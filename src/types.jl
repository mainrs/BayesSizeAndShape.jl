
#### #### #### #### #### #### 
#### Distributions
#### #### #### #### #### #### 
struct NoPriorBeta <: ContinuousUnivariateDistribution

    function NoPriorBeta()
        new()
    end

end

struct NoPriorSigma <: ContinuousMatrixDistribution

    function NoPriorSigma()
        new()
    end

end

#### #### #### #### #### #### 
#### Sample R
#### #### #### #### #### #### 
abstract type Samplermat end
struct dosamplermat <: Samplermat

    function dosamplermat()
        new()
    end

end

struct donotsamplermat <: Samplermat

    function donotsamplermat()
        new()
    end

end


#### #### #### #### #### #### 
#### dimension
#### #### #### #### #### #### 
abstract type Valuep end
struct Valuep2 <: Valuep

    function Valuep2()
        new()
    end

end

struct Valuep3 <: Valuep

    function Valuep3()
        new()
    end

end

#### #### #### #### #### #### 
#### REflections
#### #### #### #### #### ####

abstract type Reflection end
struct KeepReflection <: Reflection

    function KeepReflection()
        new()
    end

end

struct donotKeepReflection <: Reflection

    function donotKeepReflection()
        new()
    end

end

#### #### #### #### #### #### 
#### Sigma
#### #### #### #### #### ####

abstract type SigmaType end
struct GeneralSigma <: SigmaType

    function GeneralSigma()
        new()
    end

end

