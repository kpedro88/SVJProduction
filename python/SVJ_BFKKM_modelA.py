import FWCore.ParameterSet.Config as cms

processParameters = cms.vstring(
    # from svjHelper
    'HiddenValley:ffbar2Zv = on',
    '4900023:m0 = 1000',
    '4900023:mMin = 999',
    '4900023:mMax = 1001',
    '4900023:mWidth = 0.01',
    '51:m0 = 0.0',
    '51:isResonance = false',
    '53:m0 = 0.0',
    '53:isResonance = false',
    '4900211:oneChannel = 1 1 0 51 -51',
    '4900213:oneChannel = 1 1 0 53 -53',
    # from theorist card
    '4900023:0:bRatio = 0.0263',
    '4900023:1:bRatio = 0.0263',
    '4900023:2:bRatio = 0.0263',
    '4900023:3:bRatio = 0.0263',
    '4900023:4:bRatio = 0.0263',
    '4900023:5:bRatio = 0.0262',
    '4900023:6:bRatio = 0.0',
    '4900023:7:bRatio = 0.0',
    '4900023:8:bRatio = 0.0',
    '4900023:9:bRatio = 0.0',
    '4900023:10:bRatio = 0.0',
    '4900023:11:bRatio = 0.0',
    '4900023:offIfAny = 11 12 13 14 15 16',
    '4900023:12:bRatio = 0.421',
    '4900023:addChannel = 1 0.421 102 4900102 -4900102',
    '4900111:m0 = 4.0',
    '4900211:m0 = 4.0',
    '4900113:m0 = 5.0',
    '4900213:m0 = 5.0',
    'HiddenValley:FSR = on',
    'HiddenValley:fragment = on',
    'HiddenValley:alphaOrder = 1',
    'HiddenValley:Lambda = 5.0',
    'HiddenValley:pTminFSR = 5.5',
    'HiddenValley:nFlav = 2',
    'HiddenValley:spinFv = 0',
    'HiddenValley:probVector = 0.5', # rinv = 1-probVector + 0.5*probVector
    '4900113:mWidth = 3.9510605541502955e-13',
    '4900113:0:onMode = 1',
    '4900113:1:onMode = 1',
    '4900113:2:onMode = 1',
    '4900113:3:onMode = 1',
    '4900113:4:onMode = 1',
    '4900113:5:onMode = 0',
    '4900113:6:onMode = 0',
    '4900113:7:onMode = 0',
    '4900113:8:onMode = 0',
    '4900113:9:onMode = 0',
    '4900113:10:onMode = 0',
    '4900113:0:bRatio = 0.252',
    '4900113:1:bRatio = 0.252',
    '4900113:2:bRatio = 0.252',
    '4900113:3:bRatio = 0.245',
    '4900113:4:bRatio = 0.0',
    '4900113:5:bRatio = 0.0',
    '4900113:6:bRatio = 0.0',
    '4900113:7:bRatio = 0.0',
    '4900113:8:bRatio = 0.0',
    '4900113:9:bRatio = 0.0',
    '4900113:10:bRatio = 0.0',
)
