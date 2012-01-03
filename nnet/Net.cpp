/**********************************************************************
Implementation of the class handling the entire network. The
training algorithm is mainly contained here.

For latest version see: http://moonflare.com/code/nnetwork.php

Copyright (c) 2003, Derrick Coetzee
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

- The name of Derrick Coetzee may not be used to endorse or promote
  products derived from this software without specific prior
  written permission.

This software is provided by the copyright holders and contributors
"as is" and any express or implied warranties, including, but not
limited to, the implied warranties of merchantability and fitness
for a particular purpose are disclaimed. In no event shall the
copyright owner or contributors be liable for any direct, indirect,
incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or
services; loss of use, data, or profits; or business interruption)
however caused and on any theory of liability, whether in contract,
strict liability, or tort (including negligence or otherwise)
arising in any way out of the use of this software, even if advised
of the possibility of such damage.
**********************************************************************/


#include "Net.h"
#include <stdio.h>

#ifdef NEURAL_NET_DEBUG
#include <iostream>
#endif

using namespace NeuralNetwork;

// Creates a new feed-forward neural network with the given number of
// layers with the specified number of nodes in each, and the learning
// rate, momentum factor, and gain of the sigmoid function.
Net::Net(int layers, int layerSizes[],
         real initLearningRate = 0.25,
         real initMomentumFactor = 0.9,
         real initGain = 1.0,
         int max_cycles= 1000)
    :momentumFactor(initMomentumFactor),
     learningRate(initLearningRate),
     gain(initGain),max_plain_cycles(max_cycles)
{
    // Allocate and initialise layers

    numLayers = layers;
    layer = new NetLayer*[layers];
    layer[0] = new NetLayer(layerSizes[0], 0);
    for (int i=1; i<layers; i++)
    {
        layer[i] = new NetLayer(layerSizes[i], layer[i-1]);
    }

    inputLayer  = layer[0];
    outputLayer = layer[layers-1];

    // Allocate these here to avoid massive allocation slowdowns later
    int inputSize = inputLayer->getUnits();
    input = new real[inputSize];
    int outputSize = outputLayer->getUnits();
    actualOutput   = new real[outputSize];
    expectedOutput = new real[outputSize];

    // Let weights initially be random
    randomizeWeights();
    sensitivity.resize(layer[0]->units);
    stop = false;
}

// Loads a network from a stream where it was previous saved with save()
Net::Net(std::istream& in)
{
    // Allocate and initialise layers
    readRaw(in, numLayers);
    readRaw(in, momentumFactor);
    readRaw(in, learningRate);
    readRaw(in, gain);

    layer = new NetLayer*[numLayers];
    layer[0] = new NetLayer(in, 0);
    for (int i=1; i<numLayers; i++)
    {
        layer[i] = new NetLayer(in, layer[i-1]);
    }

    inputLayer  = layer[0];
    outputLayer = layer[numLayers-1];

    // Allocate these here to avoid massive allocation slowdowns later
    int inputSize = inputLayer->getUnits();
    input = new real[inputSize];
    int outputSize = outputLayer->getUnits();
    actualOutput   = new real[outputSize];
    expectedOutput = new real[outputSize];
}
// Loads a network from a stream where it was previous saved with save2()
Net::Net(FILE *fp)
{
    // Allocate and initialise layers

    fscanf(fp,"\nNumber of layers:%d\n", &numLayers);
    fscanf(fp,"Momentum factor:%lf\n",&momentumFactor);
    fscanf(fp,"Learning rate:%lf\n",&learningRate);
    fscanf(fp,"Gain:%lf\n",&gain);

    layer = new NetLayer*[numLayers];
    layer[0] = new NetLayer(fp, 0);
    for (int i=1; i<numLayers; i++)
    {
        layer[i] = new NetLayer(fp, layer[i-1]);
    }

    inputLayer  = layer[0];
    outputLayer = layer[numLayers-1];

    // Allocate these here to avoid massive allocation slowdowns later
    int inputSize = inputLayer->getUnits();
    input = new real[inputSize];
    int outputSize = outputLayer->getUnits();
    actualOutput   = new real[outputSize];
    expectedOutput = new real[outputSize];
}


// Saves network to a stream, to be read back in with the istream& constructor
void Net::save(std::ostream& out)
{
    writeRaw(out, numLayers);
    writeRaw(out, momentumFactor);
    writeRaw(out, learningRate);
    writeRaw(out, gain);

    // Save each layer
    for (int i=0; i<numLayers; i++)
    {
        layer[i]->save(out);
    }
}

// Saves network to a stream, to be read back in with the human's eye
void Net::save2()
{
    std::cout << "Number of layers: " << numLayers << std::endl;
    std::cout << "Momentum factor: " << momentumFactor << std::endl;
    std::cout << "Learning rate: " << learningRate << std::endl;
    std::cout << "Gain: " << gain << std::endl;

    // Save each layer
    for (int i=0; i<numLayers; i++)
    {
        layer[i]->save2();
    }
}




// Frees all memory allocated for network
Net::~Net()
{
    for (int i=0; i<numLayers; i++)
    {
        delete layer[i];
    }
    delete[] layer;

    delete[] input;
    delete[] actualOutput;
    delete[] expectedOutput;
}

// Initializes all weights in network to random values
void Net::randomizeWeights()
{
    // Simply call randomizeWeights on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->randomizeWeights();
    }
}

// Initializes all weights in network to zero
void Net::clearWeights()
{
    // Simply call randomizeWeights on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->clearWeights();
    }
}


// Saves weights for later restoring by restoreWeights. Only one
// set of weights can be saved, usually the best seen so far.
void Net::saveWeights()
{
    // Simply call saveWeights on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->saveWeights();
    }
}

// Restores the set of weights most recently saved by saveWeights.
void Net::restoreWeights()
{
    // Simply call restoreWeights on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->restoreWeights();
    }
}

// Propagates inputs of the net all the way through to outputs
void Net::propagate()
{
    // Simply call propagate on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->propagate(gain);
    }
}

void Net::propagate_prognosis()
{
    // Simply call propagate on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->propagate_prognosis(gain);
    }
}


// Backpropagates output errors all the way back through the network
void Net::backpropagate()
{
    // Simply call backpropagate on the layers in reverse order,
    // except the first, thus driving error back from output to input.
    for (int layerNum=numLayers-1; layerNum > 0; layerNum--) {
        layer[layerNum]->backpropagate(gain);
    }
}

void Net::_backpropagate_sens()
{
    // Simply call backpropagate on the layers in reverse order,
    // except the first, thus driving error back from output to input.
    // for sensitivity with the first!
    for (int layerNum=numLayers-1; layerNum > 0; layerNum--) {
        layer[layerNum]->backpropagate(gain);
    }

    for(int i=0; i<layer[0]->units; ++i)
        sensitivity.at(i) += layer[0]->error[i];
}

// Computes and stores error of output layer and each of its components
void Net::computeOutputError(real* target)
{
    // Just ask outputLayer to compute its error against target
    error = outputLayer->computeError(gain, target);
}

// Adjusts all weights in network to decrease error recorded by
// previous calls to backpropagate or computeOutputError.
void Net::adjustWeights()
{
    // Simply call adjustWeights on each layer except the first
    for (int layerNum=1; layerNum < numLayers; layerNum++) {
        layer[layerNum]->adjustWeights(momentumFactor, learningRate);
    }
}

// Sets the values of the inputs to the network
void Net::setInputs(real* inputs)
{
    // Sets output values of input-layer perceptrons, which is input
    // of network
    inputLayer->setOutputs(inputs);
}

// Gets the outputs and places them in the array outputs
void Net::getOutputs(real* outputs)
{
    // Gets output of output-layer perceptrons, which is output of
    // network
    outputLayer->getOutputs(outputs);
}

// Trains a single training example once
void Net::simpleTrain(real* input, real* expectedOutput)
{
    // See first what the network produces now for the input
    // and see how far off it is.
    setInputs(input);
    propagate();
    computeOutputError(expectedOutput);

    // Backpropagate that error data and then adjust the weights based
    // on it to reduce total error as quickly as possible.
    backpropagate();
    adjustWeights();
}

// Trains a single sensitivity training example once
void Net::_simpleTrain_sens(real* input, real* expectedOutput)
{
    // See first what the network produces now for the input
    // and see how far off it is.
    setInputs(input);
    propagate();
    computeOutputError(expectedOutput);

    // Backpropagate that error data and then adjust the weights based
    // on it to reduce total error as quickly as possible.

    _backpropagate_sens();

}


// Trains for a given number of epochs on an entire training set
void Net::train(int epochs, ExampleFactory &trainingExamples)
{
    int inputSize = inputLayer->getUnits();
    int outputSize = outputLayer->getUnits();

    // Train on each training example an average of epochs times
    for (int n=0; n < epochs*trainingExamples.numExamples(); n++) {
        trainingExamples.getExample(inputSize, input, outputSize, expectedOutput);
        simpleTrain(input, expectedOutput);
    }
}

// Trains for a given number of epochs on an entire training set
void Net::_train_sens(int epochs, ExampleFactory &trainingExamples)
{

    int inputSize = inputLayer->getUnits();
    int outputSize = outputLayer->getUnits();

    // Train on each training example an average of epochs times
    for (int n=0; n < epochs*trainingExamples.numExamples(); n++) {
        trainingExamples.getExample(inputSize, input, outputSize, expectedOutput);
        _simpleTrain_sens(input, expectedOutput);
    }
}


// Tests the network using the given examples, returning total error
real Net::test(ExampleFactory &testExamples)
{
    int inputSize = inputLayer->getUnits();
    int outputSize = outputLayer->getUnits();

    real totalError = 0;

    // Run network once on each example, adding error each time to a
    // running total.
    for (int n=0; n < testExamples.numExamples(); n++) {
        testExamples.getExample(inputSize, input, outputSize, expectedOutput);
        run(input, actualOutput);
        computeOutputError(expectedOutput);
        totalError += error;
    }

#ifdef NEURAL_NET_DEBUG
//    std::cout << "Error: " << totalError << std::endl;
#endif

    return totalError;
}

// Automatically trains the network until its performance on
// the test set appears to achieve a maximum. Returns total
// error on the test set after completion.
// - epochsBetweenTests establishes how many tests are done between
//   test set checks for accuracy.
// - cutOffError establishes how much worse, as a multiple, error has
//   to be than the minimum error seen before we stop training. Must be >1.
real Net::autotrain(ExampleFactory &trainingExamples,
                    ExampleFactory &testExamples,
                    int epochsBetweenTests,
                    float cutOffError)
{
    // Get initial error with current weight set
    real minTestError = test(testExamples);
    real testError = minTestError;

    int cycle=0;

    while (cycle < max_plain_cycles && !stop) {
        cycle++;
//	std::cout << cycle << std::endl;
        // Train for a while on the training examples
        train(epochsBetweenTests, trainingExamples);

        // How good is network now? Save weights if it's the best
        // we've seen so far on the test set.

        testError = test(testExamples);
        std::cout << testError << " " << cycle << std::endl;

        if (testError < minTestError) {
            saveWeights();
            minTestError = testError;
            cycle = 0;
        }
    }

    // Restore weights so performance on test set is best we ever saw
    restoreWeights();

    return minTestError;
}

real Net::autotrain_occ(ExampleFactory &trainingExamples,
                        int epochsBetweenTests,
                        float cutOffError)
{
    // Get initial error with current weight set
    real minTestError = 1.0e30;
    real testError = minTestError;

    int cycle=0;

    while (cycle < max_plain_cycles && !stop) {
        cycle++;
//	std::cout << cycle << std::endl;
        // Train for a while on the training examples
        train(epochsBetweenTests, trainingExamples);

        // How good is network now? Save weights if it's the best
        // we've seen so far on the test set.

        testError = 1.0-30;
//	std::cout << testError << " " << cycle << std::endl;

        //std::cout << testError << std::endl;

        if (testError < minTestError) {
            saveWeights();
            minTestError = testError;
            cycle = 0;
        }
    }

    // Restore weights so performance on test set is best we ever saw
    restoreWeights();

    return minTestError;
}


// Runs the network on an input and feeds it forward to produce an
// output. For usage after training with autotrain.
void Net::run(real* input, real* output)
{
    setInputs(input);
    propagate();
    getOutputs(output);
}

void Net::run_prognosis(real* input, real* output)
{
    setInputs(input);
    propagate_prognosis();
    getOutputs(output);
}



// Deallocates storage used only during training
void Net::doneTraining()
{
    // Simply call doneTraining() on each layer
    for (int layerNum=0; layerNum < numLayers; layerNum++) {
        layer[layerNum]->doneTraining();
    }
}
