local class = require 'ext.class'
local ForwardEuler = require 'hydro.int.fe'
local IterativeCrankNicolson = require 'hydro.int.icn'
local BackwardEuler = require 'hydro.int.be'
local BackwardEulerCPU = require 'hydro.int.be-cpu'
local RungeKutta = require 'hydro.int.rk'

--the following are from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method

local RungeKutta2 = class(RungeKutta)
RungeKutta2.name = 'Runge-Kutta 2'
RungeKutta2.alphas = {
	{1,0},
	{1,0},
}
RungeKutta2.betas = {
	{.5, 0},
	{0, 1},
}

local RungeKutta2Heun = class(RungeKutta)
RungeKutta2Heun.name = 'Runge-Kutta 2 Heun'
RungeKutta2Heun.alphas = { 
	{1, 0},
	{1, 0},
}
RungeKutta2Heun.betas = {
	{1, 0},
	{.5, .5},
}

local RungeKutta2Ralston = class(RungeKutta)
RungeKutta2Ralston.name = 'Runge-Kutta 2 Ralston'
RungeKutta2Ralston.alphas = {
	{1, 0},
	{1, 0},
}
RungeKutta2Ralston.betas = {
	{2./3., 0},
	{1./4., 3./4.},
}

local RungeKutta3 = class(RungeKutta)
RungeKutta3.name = 'Runge-Kutta 3'
RungeKutta3.alphas = {
	{1, 0, 0},
	{1, 0, 0},
	{1, 0, 0},
}
RungeKutta3.betas = {
	{.5, 0, 0},
	{-1, 2, 0},
	{1./6., 2./6., 1./6.},
}

local RungeKutta4 = class(RungeKutta)
RungeKutta4.name = 'Runge-Kutta 4'
RungeKutta4.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
}
RungeKutta4.betas = {
	{.5, 0, 0, 0},
	{0, .5, 0, 0},
	{0, 0, 1, 0},
	{1./6., 2./6., 2./6., 1./6.},
}

local RungeKutta4_3_8thsRule = class(RungeKutta)
RungeKutta4_3_8thsRule.name = 'Runge-Kutta 4, 3/8ths rule'
RungeKutta4_3_8thsRule.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
}
RungeKutta4_3_8thsRule.betas = {
	{1./3., 0, 0, 0},
	{-1./3., 0, 0, 0},
	{1, -1, 1, 0},
	{1./8., 3./8., 3./8., 1./8.},
}

--the following are from 1998 Gottleib, Shu "Total Variation Diminishing Runge-Kutta Methods"

local RungeKutta2TVD = class(RungeKutta)
RungeKutta2TVD.name = 'Runge-Kutta 2, TVD'
RungeKutta2TVD.alphas = {
	{1, 0},
	{.5, .5},
}
RungeKutta2TVD.betas = {
	{1, 0},
	{0, .5},
}

local RungeKutta2NonTVD = class(RungeKutta)
RungeKutta2NonTVD.name = 'Runge-Kutta 2, non-TVD'
RungeKutta2NonTVD.alphas = {
	{1, 0},
	{1, 0},
}
RungeKutta2NonTVD.betas = {
	{-20, 0},
	{41./40., -1./40.},
}

local RungeKutta3TVD = class(RungeKutta)
RungeKutta3TVD.name = 'Runge-Kutta 3, TVD'
RungeKutta3TVD.alphas = {
	{1, 0, 0},
	{3/4, 1/4, 0},
	{1/3, 0, 2/3},
}
RungeKutta3TVD.betas = {
	{1, 0, 0},
	{0, 1/4, 0},
	{0, 0, 2/3},
}

local RungeKutta4TVD = class(RungeKutta)
RungeKutta4TVD.name = 'Runge-Kutta 4, TVD'
RungeKutta4TVD.alphas = {
	{1, 0, 0, 0},
	{649./1600., 951./1600., 0, 0},
	{53989./2500000., 4806213./20000000., 23619./32000., 0},
	{1./5., 6127./30000., 7873./30000., 1./3.},
}
RungeKutta4TVD.betas = {
	{.5, 0, 0, 0},
	{-10890423./25193600., 5000./7873, 0, 0},
	{-102261./5000000., -5121./20000., 7873./10000., 0},
	{1./10., 1./6., 0, 1./6.},
}

-- this one is from 1995 Jiang, Shu "Efficient Implementation of Weighted ENO Schemes"
--  found at http://lsec.cc.ac.cn/lcfd/DEWENO/paper/WENO_1996.pdf

local RungeKutta4NonTVD = class(RungeKutta)
RungeKutta4NonTVD.name = 'Runge-Kutta 4, non-TVD'
RungeKutta4NonTVD.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{-1./3., 1./3., 2./3., 1./3.},
}
RungeKutta4NonTVD.betas = {
	{.5, 0, 0, 0},
	{0, .5, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1./6.},
}

return require 'ext.table'{
	ForwardEuler,
	IterativeCrankNicolson,
	RungeKutta2,
	RungeKutta2Heun,
	RungeKutta2Ralston,
	RungeKutta3,
	RungeKutta4,
	RungeKutta4_3_8thsRule,
	RungeKutta2TVD,
	RungeKutta2NonTVD,
	RungeKutta3TVD,
	RungeKutta4TVD,
	RungeKutta4NonTVD,
	BackwardEuler,
	BackwardEulerCPU,
}
