#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for the workflow_with_facets.cwl using medium sized dataset

Takes about  9467.131s to complete

Usage:

    $ INTEGRATION_TESTS=True USE_LSF=True CWL_ENGINE=toil python3 tests/test_workflow_with_facets.medium.py

"""
import os
import sys
import unittest



from pluto import (
    PlutoTestCase,
    ENABLE_INTEGRATION_TESTS,
    OFile,
    ODir
)


try:
    from fixtures import WORKFLOW_MEDIUM_JSON
except ModuleNotFoundError:
    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
    from fixtures import WORKFLOW_MEDIUM_JSON
    

class TestWorkflowWithFacetsMedium(PlutoTestCase):
    cwl_file = 'workflow_with_facets.cwl'

    @unittest.skipIf(ENABLE_INTEGRATION_TESTS!=True, "is a large integration test")
    def test_medium_workflow(self):
        """
        ...
        ---------
        Ran 1 test in 8698.861s
        """
        input_json = WORKFLOW_MEDIUM_JSON
        output_json, output_dir = self.run_cwl(input = input_json, input_is_file = True)
        expected_output = {
        'analysis_dir': ODir(name='analysis', items=[
                OFile(name='11089_G.gene.cna.txt', size=6547460, hash='72c2eb26e69f8dcb53cb0ff2d90997f2ae69f308'),
                OFile(name='11089_G.muts.maf', size=416363, hash='16ee2715491a5f61810cfd7233163905fb497361'),
                OFile(name='11089_G.muts.share.maf', size=66088, hash='8b3e0924733421535a3fc5a5d174eef4d5d92a28'),
                OFile(name='11089_G.seg.cna.txt', size=69171, hash='7d13e1422b9d7112e60f9ab748d524d7446d1dde'),
                OFile(name='11089_G.svs.maf', size=706820, hash='a4a799461650e4cc832e54de2df50a9e1b28a41c')], dir=output_dir),
        'facets_dir': ODir(name='facets', items=[
            ODir(name='Sample21.s_C_Patient21_N001_d', items=[
                OFile(name='Sample21.s_C_Patient21_N001_d_hisens.ccf.portal.maf', size=108887494, hash='f4f5d6e0dae3c8ad89279910a4e083dc8471bcc2'),
                OFile(name='Sample21.arm_level.txt', size=1900, hash='8af0ec2771d47ecb9c17727755094474d32a9f8b'),
                OFile(name='Sample21.txt', size=591, hash='4350cbbdf4151b1f6a26b5269d4a227faed40f06'),
                OFile(name='Sample21.gene_level.txt', size=569488, hash='6e663d203f0b17b52021c93e1b10cdcef6bf88b6'),
                OFile(name='Sample21_hisens.cncf.txt', size=3642, hash='0ea89011191860974a6e98bc0fc8a55752979856'),
                OFile(name='Sample21_hisens.rds', size=616659, hash='696afa887ffe3ea56de3824ab2d9483a8878186a'),
                OFile(name='Sample21_hisens.seg', size=1361, hash='e62dd226d7d14884146ca263a7bf689108ebea7b'),
                OFile(name='Sample21_purity.rds', size=616544, hash='bd51898f1d0d04bd761cc9e7ada7533c04f8cc7c'),
                OFile(name='Sample21_purity.seg', size=1258, hash='b4fb384eb9dd63a938fe8b38001ad6b4a38fc8c5'),
                OFile(name='Sample21.qc.txt', size=1283, hash='cfe0031e6ca74c700f7970839252376dda907deb')]),
            ODir(name='Sample5.s_C_Patient5_N001_d', items=[
                OFile(name='Sample5.s_C_Patient5_N001_d_hisens.ccf.portal.maf', size=51378563, hash='1886433b354d192591424faba968c975e880421b'),
                OFile(name='Sample5.arm_level.txt', size=1504, hash='ce1d534e85fe61b72d7d732f0c5b94b9be4cd5bc'),
                OFile(name='Sample5.txt', size=594, hash='df5143e3ef42c3546dd3a53ddd7c69679c94ece8'),
                OFile(name='Sample5.gene_level.txt', size=701138, hash='107f1665e94907ca8ac507852b7031458a1524af'),
                OFile(name='Sample5_hisens.cncf.txt', size=11955, hash='0e3745e5c13545ac9e179fc18b093f2197f095f0'),
                OFile(name='Sample5_hisens.rds', size=768703, hash='1508ccc6b3334e14dcc25439668689fd8416c590'),
                OFile(name='Sample5_hisens.seg', size=4206, hash='81d558ac13df8f87edb02759b5aa01a6a8065904'),
                OFile(name='Sample5_purity.rds', size=766511, hash='5adc9fb028bcb7752a7550cc3c2e180a032179c7'),
                OFile(name='Sample5_purity.seg', size=2285, hash='a1a35082d279dbc8cb06ad1a6f0b42a1d5601f1e'),
                OFile(name='Sample5.qc.txt', size=1386, hash='540eb51045962627fe6c17e47b0f69e3cb35e081')]),
            ODir(name='Sample1.s_C_Patient1_N001_d', items=[
                OFile(name='Sample1.s_C_Patient1_N001_d_hisens.ccf.portal.maf', size=174888292, hash='0f9f3acc5b9bfeedc35713635eba030adf15e690'),
                OFile(name='Sample1.arm_level.txt', size=1629, hash='07f5b70585370c5926f2849995f5d74af223e0db'),
                OFile(name='Sample1.txt', size=590, hash='02018eabc346542db0589162f9a3a073726a3dd6'),
                OFile(name='Sample1.gene_level.txt', size=901681, hash='a62c291ff0bc23c12cb1939a545925b3192ed2ea'),
                OFile(name='Sample1_hisens.cncf.txt', size=17029, hash='0dc6f9ccb68f24ee50ff7cbe26ec37b537898625'),
                OFile(name='Sample1_hisens.rds', size=728644, hash='5f64eb5770316eba7d1df10db097b2ee7ba609d9'),
                OFile(name='Sample1_hisens.seg', size=5691, hash='64d2d15fc606b52f78d724eaa16aef05783bbf32'),
                OFile(name='Sample1_purity.rds', size=725530, hash='26acecdeaa02189b43cb4ef8f747dca04d97f21b'),
                OFile(name='Sample1_purity.seg', size=2967, hash='5f28d741b409e7dc92f9a2b4d502e657c56f8f9e'),
                OFile(name='Sample1.qc.txt', size=1361, hash='fd3593a93588eb43ed78b74978704a5de3e56768')]),
            ODir(name='Sample12.s_C_Patient12_N001_d',items=[
                OFile(name='Sample12.s_C_Patient12_N001_d_hisens.ccf.portal.maf', size=33858338, hash='1bc7c99356930074ee412383a8798fec0ddea103'),
                OFile(name='Sample12.arm_level.txt', size=1515, hash='47f866d25f6f2b2d226b64c76fc6b6d08ea46658'),
                OFile(name='Sample12.txt', size=762, hash='07d1ea16f4ba78e59afc73c519cc95235053c502'),
                OFile(name='Sample12.gene_level.txt', size=939824, hash='b6b7b1b792519139ad889b654e0fc4ac96486ae0'),
                OFile(name='Sample12_hisens.cncf.txt', size=14938, hash='5e12530c6e9baebfbcaccecf55bd3cd34efd73ba'),
                OFile(name='Sample12_hisens.rds', size=905882, hash='5b188c566a89510c3c5d2901fec48973226be27e'),
                OFile(name='Sample12_hisens.seg', size=4805, hash='ead3884d2ac67dc7f58299d69f3c79df2a34b942'),
                OFile(name='Sample12_purity.rds', size=903732, hash='7f1b31c5c4bb33e08b20d293852b7b90232c213f'),
                OFile(name='Sample12_purity.seg', size=2697, hash='fe306a606d349eba7116d3d9c8bc053f765c40f0'),
                OFile(name='Sample12.qc.txt', size=1368, hash='225589806bd1418b69f6e6b126f7854bc587c93f')]),
            ODir(name='Sample17.s_C_Patient17_N001_d', items=[
                OFile(name='Sample17.s_C_Patient17_N001_d_hisens.ccf.portal.maf', size=61845296, hash='8cdd2bc33f12565429c4d550dd23b871f628b380'),
                OFile(name='Sample17.arm_level.txt', size=1451, hash='b040ff1705a708b9c3532f68af9b7ddce9c2e1af'),
                OFile(name='Sample17.txt', size=711, hash='30c3a6df33fda7a4660e479ec83922b00a6f36dc'),
                OFile(name='Sample17.gene_level.txt', size=674625, hash='36f45049e0a4f09e84477fb013cce104229fbe67'),
                OFile(name='Sample17_hisens.cncf.txt', size=13569, hash='3a9e2d666c3dcd9ec9f6958f9e55b53a028e206b'),
                OFile(name='Sample17_hisens.rds', size=718061, hash='621e576c38bb376b79e46fa2a327aede80933e5d'),
                OFile(name='Sample17_hisens.seg', size=4397, hash='b085804e5608e5577623ec8fe9d233b781d3b644'),
                OFile(name='Sample17_purity.rds', size=716807, hash='dbff46198bea8c4cdad0f1a65391e94eb37be0ff'),
                OFile(name='Sample17_purity.seg', size=3066, hash='eff585106b1ff3e857ed074c241885af3844859b'),
                OFile(name='Sample17.qc.txt', size=1353, hash='8bde8eb16eb3eb59d2a8b03d5f129ab649149f65')]),
            ODir(name='Sample28.s_C_Patient28_N001_d', items=[
                OFile(name='Sample28.s_C_Patient28_N001_d_hisens.ccf.portal.maf', size=120654677, hash='0c4ac1ca2c450fc16213a9aa1515b594c24814d0'),
                OFile(name='Sample28.arm_level.txt', size=1909, hash='f9f9f8880538fd7833781f0a472d06c2e68cfc6c'),
                OFile(name='Sample28.txt', size=597, hash='f553ef7868018ff066eae92fec5db08a19a4f197'),
                OFile(name='Sample28.gene_level.txt', size=509340, hash='8d198a5658804849c9d99b210a6e0615eda23a67'),
                OFile(name='Sample28_hisens.cncf.txt', size=4623, hash='8b4e4339705616a2eb9a611d26668ff259e0ea93'),
                OFile(name='Sample28_hisens.rds', size=581539, hash='71b21ab0abf91eb2a31a19b1c36ca4ff4fe0ab82'),
                OFile(name='Sample28_hisens.seg', size=1712, hash='870a4e66276768736191a622523b5664dfd5823b'),
                OFile(name='Sample28_purity.rds', size=581286, hash='06d334d5ab287d12aef798cfad5d87f56084ce54'),
                OFile(name='Sample28_purity.seg', size=1413, hash='1f4df7ff87caa5b2c39922d795a109f4374b5a34'),
                OFile(name='Sample28.qc.txt', size=1329, hash='c9f96b6e74ad2aad0aaebef8a2f9fb0c1df8cbd7')]),
            ODir(name='Sample16.s_C_Patient16_N001_d', items=[
                OFile(name='Sample16.s_C_Patient16_N001_d_hisens.ccf.portal.maf', size=46188384, hash='542851d8d78ae2022d2d319298f2fca7fa1d94cb'),
                OFile(name='Sample16.arm_level.txt', size=1709, hash='e5bd70471da5aceaf545be6405c5ab44da91064f'),
                OFile(name='Sample16.txt', size=601, hash='8002c18881cecff90e3762d013974fb6981ea75b'),
                OFile(name='Sample16.gene_level.txt', size=600614, hash='d4881c3f65a84e73195e8b56de325dfcf4e4e01a'),
                OFile(name='Sample16_hisens.cncf.txt', size=8976, hash='147e638a415fc1c904d02431b2874efe33d0b342'),
                OFile(name='Sample16_hisens.rds', size=773030, hash='59bab59e78318eb267d9589aa6f650959b47bc6a'),
                OFile(name='Sample16_hisens.seg', size=3334, hash='1bb1490d894699925fc1f5ab89e4cbcb568f82b9'),
                OFile(name='Sample16_purity.rds', size=771948, hash='f1b173c8507ad3b2a5dbf4cc9a95902bf408abd4'),
                OFile(name='Sample16_purity.seg', size=2191, hash='b4065c5603c21724e5e74c8d49e2a0c97f712f2c'),
                OFile(name='Sample16.qc.txt', size=1358, hash='5c982742449b705c0224c0656fcfea17baf46ee1')]),
            ODir(name='Sample11.s_C_Patient11_N001_d', items=[
                OFile(name='Sample11.s_C_Patient11_N001_d_hisens.ccf.portal.maf', size=36352097, hash='675e67d7ddc6ed75a41c458481bee8aa5bae857d'),
                OFile(name='Sample11.arm_level.txt', size=1525, hash='13a74ad441b2c24a4affb855f8a79483846bfad5'),
                OFile(name='Sample11.txt', size=639, hash='20a56017fd43021a02312fba1762420f40dc6c85'),
                OFile(name='Sample11.gene_level.txt', size=671056, hash='bff64a96cf49f180086edf916e6c6bdf603b5ef3'),
                OFile(name='Sample11_hisens.cncf.txt', size=10363, hash='3827502368c7347a3924a716b6fe60544cb3cc2f'),
                OFile(name='Sample11_hisens.rds', size=720072, hash='bbd4d3725cee4707c957ef4aff8ea5ef811c189f'),
                OFile(name='Sample11_hisens.seg', size=3638, hash='75f2c0c83b1ac1cc8cc9bf819c71ac504fa9bee2'),
                OFile(name='Sample11_purity.rds', size=719069, hash='17bcd89e06710f4086f87f023c6a182e493e3cc3'),
                OFile(name='Sample11_purity.seg', size=2744, hash='ce3a31c4906a962ee2b43976132ef538168a697e'),
                OFile(name='Sample11.qc.txt', size=1352, hash='43f0294dfb40a0421e56ea35926bfda484d43d23')]),
            ODir(name='Sample13.s_C_Patient13_N001_d', items=[
                OFile(name='Sample13.s_C_Patient13_N001_d_hisens.ccf.portal.maf', size=27544622, hash='f9785eee1574f6e85df67b80f6b60879fb7a1a61'),
                OFile(name='Sample13.arm_level.txt', size=1619, hash='2c1211c420b4da7f336bd2fe99156e3b9b0e7c32'),
                OFile(name='Sample13.txt', size=672, hash='b182fc5f96d497cdfced581afd409ab78e34eba8'),
                OFile(name='Sample13.gene_level.txt', size=882343, hash='c8f71dfcd5244f344c1a481864f924eeb84ede65'),
                OFile(name='Sample13_hisens.cncf.txt', size=9876, hash='3a9b5c3ada3eb41562781093ccdd211f43ddee96'),
                OFile(name='Sample13_hisens.rds', size=923512, hash='f2db4d71b1d9a1b8da90603a81bbf514337869e2'),
                OFile(name='Sample13_hisens.seg', size=3501, hash='05aef3b06f91eac0bc8a2034ffa75a2b585c5eeb'),
                OFile(name='Sample13_purity.rds', size=921728, hash='3bf8506d3d77ad2a2cce4e39bbbbaed2bb5e3343'),
                OFile(name='Sample13_purity.seg', size=1872, hash='27e67ba4902078d903b2317933bbbb0e36f16736'),
                OFile(name='Sample13.qc.txt', size=1351, hash='71f68324b72f77e5d494c51791e08bccb6541467')]),
            ODir(name='Sample3.s_C_Patient3_N001_d', items=[
                OFile(name='Sample3.s_C_Patient3_N001_d_hisens.ccf.portal.maf', size=56270406, hash='d0115fa74a43e57498c11eb5148d8e7ef3c03fc9'),
                OFile(name='Sample3.arm_level.txt', size=1865, hash='6ebe1352c2b63060f3dfc6ceb97a499dc04f318a'),
                OFile(name='Sample3.txt', size=578, hash='270230a33f6571484d525b11214b1f2471b2dd56'),
                OFile(name='Sample3.gene_level.txt', size=706751, hash='87320ad1d73e0a1e0ea52cf83ca1c24a18d94f40'),
                OFile(name='Sample3_hisens.cncf.txt', size=3343, hash='38294c690440e183e8b4488a4452edd0ce60fab2'),
                OFile(name='Sample3_hisens.rds', size=815782, hash='ec3e53d6131de9119ea20b59b32919574dff1679'),
                OFile(name='Sample3_hisens.seg', size=1237, hash='8816010bda1aa473b55f972e5bf9302ab8877c87'),
                OFile(name='Sample3_purity.rds', size=815832, hash='b47f2833cb6f00e5c144fa28249f69434e9a76ea'),
                OFile(name='Sample3_purity.seg', size=1237, hash='8816010bda1aa473b55f972e5bf9302ab8877c87'),
                OFile(name='Sample3.qc.txt', size=1254, hash='4981ecfa24438f1b7d622e745ba5020eaf60fa78')]),
            ODir(name='Sample7.s_C_Patient7_N001_d', items=[
                OFile(name='Sample7.s_C_Patient7_N001_d_hisens.ccf.portal.maf', size=45465885, hash='d6f82306ee11ec48d54656119ed1272d3ef25e3b'),
                OFile(name='Sample7.arm_level.txt', size=1840, hash='c8a7436d18e42ed23011324ddfbe2df083b3ca79'),
                OFile(name='Sample7.txt', size=710, hash='988de61c55ab77422f29d1341b2ef865dbea5f60'),
                OFile(name='Sample7.gene_level.txt', size=553685, hash='bd0caff5b6612bbfc8bb8771965ccf9b631c8675'),
                OFile(name='Sample7_hisens.cncf.txt', size=7995, hash='c520b7908d699e27f05e6e3c62dffa54aad369ff'),
                OFile(name='Sample7_hisens.rds', size=622782, hash='038d4df531117a3f9d0113925d9b69df9b190c9c'),
                OFile(name='Sample7_hisens.seg', size=2674, hash='6e7dbb031f046952400f0645b0fe93476e580c74'),
                OFile(name='Sample7_purity.rds', size=622161, hash='05e660d9c69275f0984690b9e327782d3a92b7b9'),
                OFile(name='Sample7_purity.seg', size=1913, hash='5df4463fe377d46e0b4d1714bae9a152ec094c86'),
                OFile(name='Sample7.qc.txt', size=1337, hash='5b5b331eee7a5f05a1c8b997949087fba214d161')]),
            ODir(name='Sample10.s_C_Patient10_N001_d', items=[
                OFile(name='Sample10.s_C_Patient10_N001_d_hisens.ccf.portal.maf', size=63209220, hash='3600ec94d0321619b3294437b91a1b8beb79c97b'),
                OFile(name='Sample10.arm_level.txt', size=1887, hash='637fa7259d0fa7e0d4c43ef6c68209c63c62b2dd'),
                OFile(name='Sample10.txt', size=663, hash='b9ac85144c27a8d61d4042a20942d0df43152bb6'),
                OFile(name='Sample10.gene_level.txt', size=1026175, hash='692a365e399f8c5636caeed38a8dcd9b1010905b'),
                OFile(name='Sample10_hisens.cncf.txt', size=13605, hash='06520f90f7f38440d0e66bfe08874cbc97f2d0c9'),
                OFile(name='Sample10_hisens.rds', size=891327, hash='f004ffd65e33e55bb0505dc898a32189d16c24ad'),
                OFile(name='Sample10_hisens.seg', size=4591, hash='e6f9f208f8bae0a3080d221be79d422660eb0167'),
                OFile(name='Sample10_purity.rds', size=889163, hash='35f953edbdca56e089874acf4217321b7cd94944'),
                OFile(name='Sample10_purity.seg', size=2214, hash='b1c99eb39bdcd3389291e0c039fcb45aefdd7596'),
                OFile(name='Sample10.qc.txt', size=1388, hash='4127ac753a99dce3c9081021813e5784f713d8e3')]),
            ODir(name='Sample18.s_C_Patient18_N001_d', items=[
                OFile(name='Sample18.s_C_Patient18_N001_d_hisens.ccf.portal.maf', size=96929702, hash='af2afa576708903f6d9829bbcb23cbf66c73a09e'),
                OFile(name='Sample18.arm_level.txt', size=1723, hash='466b88b52d65af0dd0dcd19dda3822ca68d38d9f'),
                OFile(name='Sample18.txt', size=720, hash='df5d84236277931cdeb9543f6d37edeada3a5f0d'),
                OFile(name='Sample18.gene_level.txt', size=621856, hash='f17adc4b65d8f7b08e71465c6d1602844a4527b8'),
                OFile(name='Sample18_hisens.cncf.txt', size=8098, hash='639f929dbd8092af2fa97718877bc2a7438e0429'),
                OFile(name='Sample18_hisens.rds', size=622428, hash='c332f5db8964a9a4df66f972e2899417c1b24eee'),
                OFile(name='Sample18_hisens.seg', size=2729, hash='9bb5248ac62117eb8bfe5c6f34ae4103a844c61e'),
                OFile(name='Sample18_purity.rds', size=621691, hash='e81821a181b199c4e61185808592b36d89e79bd9'),
                OFile(name='Sample18_purity.seg', size=1895, hash='d9eac576e97b6149f964e2bd96585b0be3bf3f3c'),
                OFile(name='Sample18.qc.txt', size=1367, hash='edf6c6d86850b62951b916ca1223a43a9e641c70')]),
            ODir(name='Sample24.s_C_Patient24_N001_d', items=[
                OFile(name='Sample24.s_C_Patient24_N001_d_hisens.ccf.portal.maf', size=95307035, hash='5f8bea344bc62637449897caa4c3f5a51c51b880'),
                OFile(name='Sample24.arm_level.txt', size=1906, hash='f7963d213de2554283673bbaad2f6758768d160a'),
                OFile(name='Sample24.txt', size=584, hash='d7a9d3255cebef2d35f2e8f41ba1b1815a35af37'),
                OFile(name='Sample24.gene_level.txt', size=707804, hash='77266581d09b9b405a4d34f0540e785dbeb379fd'),
                OFile(name='Sample24_hisens.cncf.txt', size=4922, hash='a64b18a3245f569262b53f6e15aad71d67501be8'),
                OFile(name='Sample24_hisens.rds', size=786110, hash='1cc8dde029f75be5965288ea53ce44291e39777e'),
                OFile(name='Sample24_hisens.seg', size=1884, hash='073b0a35d2ce628703d1526ca5c7bac21e693a87'),
                OFile(name='Sample24_purity.rds', size=785433, hash='0af6d15d483c915cc09da4fd3dc7be129434a949'),
                OFile(name='Sample24_purity.seg', size=1257, hash='ace13a09fe8f059592449d4b0648e53623de9c9a'),
                OFile(name='Sample24.qc.txt', size=1266, hash='eab37d701dfcaf23d798bc5dd0881818ca1af4c4')]),
            ODir(name='Sample19.s_C_Patient19_N001_d', items=[
                OFile(name='Sample19.s_C_Patient19_N001_d_hisens.ccf.portal.maf', size=80099004, hash='139d09c84ce8ceb6e2a698772266d8570d671fa5'),
                OFile(name='Sample19.arm_level.txt', size=1765, hash='c8a3b83caf7fc583860b4fe241838b0c7c13f1d8'),
                OFile(name='Sample19.txt', size=600, hash='1a37025c41d74df6075ddb09bb771ce8ecca9d32'),
                OFile(name='Sample19.gene_level.txt', size=625043, hash='a2fc4a6258a5281cd99c43de2c94a65340ebc382'),
                OFile(name='Sample19_hisens.cncf.txt', size=7064, hash='a7bde8b4f196494cf86ddf501bcba7775f24e8cb'),
                OFile(name='Sample19_hisens.rds', size=669558, hash='d3a78342b1d58c716920f345a6a721032dabda1b'),
                OFile(name='Sample19_hisens.seg', size=2561, hash='3ae876d1a5ed5620b72151fdcdbe392f3c90aefd'),
                OFile(name='Sample19_purity.rds', size=668833, hash='106d18871fcf604a8c8bfa68557e2ca6585edf44'),
                OFile(name='Sample19_purity.seg', size=1880, hash='95e4b5ee4fd0c0c25dc6ed57c50e3702d27e1c00'),
                OFile(name='Sample19.qc.txt', size=1352, hash='4f9987fcea002d3867c39cdd9533a8bcf0f17c45')]),
            ODir(name='Sample20.s_C_Patient20_N001_d', items=[
                OFile(name='Sample20.s_C_Patient20_N001_d_hisens.ccf.portal.maf', size=144041600, hash='96f18c4fabdfee52beed70713caa7441e477e7b6'),
                OFile(name='Sample20.arm_level.txt', size=1912, hash='802bfddeeb394f0c896a9ed58caf65731ce60971'),
                OFile(name='Sample20.txt', size=667, hash='8ae733e7b70ec4a25568962c3f15ea0573df64dd'),
                OFile(name='Sample20.gene_level.txt', size=788731, hash='998e9f40e59cd9d77448c21de3013265348bf314'),
                OFile(name='Sample20_hisens.cncf.txt', size=7826, hash='d8c606880116aab31980d3efbfb33ef955c2acf4'),
                OFile(name='Sample20_hisens.rds', size=900050, hash='2123ac1806edd2772d669630d48ca980cf974ba8'),
                OFile(name='Sample20_hisens.seg', size=2770, hash='26d8da55158fcea4607bbc8125071921de1513f4'),
                OFile(name='Sample20_purity.rds', size=898724, hash='ca2c03dcd69f89de30c0cc34ec42208fbf075e75'),
                OFile(name='Sample20_purity.seg', size=1563, hash='6ca959b18963fc7a22300288e6caf0ea85a4b732'),
                OFile(name='Sample20.qc.txt', size=1363, hash='514e9cb3033abd02b1fbf8c5a0515e15052af6b9')]),
            ODir(name='Sample15.s_C_Patient15_N001_d', items=[
                OFile(name='Sample15.s_C_Patient15_N001_d_hisens.ccf.portal.maf', size=63778532, hash='e40e8bae617a6b1b7adc41084ec9d625d5771f9c'),
                OFile(name='Sample15.arm_level.txt', size=1824, hash='929524589e1adca5241872d9ecb0dce9897bbabe'),
                OFile(name='Sample15.txt', size=798, hash='f2fdbdb1c496322ecceeb0f6872b9e57b1d6d125'),
                OFile(name='Sample15.gene_level.txt', size=430468, hash='9d3e2c47beb125fb187b24f7e8b7afe0e738b4cb'),
                OFile(name='Sample15_hisens.cncf.txt', size=11654, hash='f5f89eb1c57f548df0fc663938258cbe7a460cba'),
                OFile(name='Sample15_hisens.rds', size=534960, hash='f6c04f62fef06207ff700181b012f57dbd713299'),
                OFile(name='Sample15_hisens.seg', size=4182, hash='1715e41998fd8dd85d1d8908073acca98660ce29'),
                OFile(name='Sample15_purity.rds', size=533755, hash='793dbd9fb47e30142aa638489357f1ea64a8273c'),
                OFile(name='Sample15_purity.seg', size=2839, hash='0a575538ee0b4a174c0639be4ed6b60212fa0a40'),
                OFile(name='Sample15.qc.txt', size=1389, hash='e99e99ec7752f4ef457a24994a2bfeb1e3d3de68')]),
            ODir(name='Sample31.s_C_Patient31_N001_d', items=[
                OFile(name='Sample31.s_C_Patient31_N001_d_hisens.ccf.portal.maf', size=42776060, hash='d2b87754e40ada818a304883b991e7580baed8b4'),
                OFile(name='Sample31.arm_level.txt', size=1859, hash='6426c07b429b4174a86adca733fab619f26e9a88'),
                OFile(name='Sample31.txt', size=597, hash='546a9a543daa2882566658c6461187dc4e36da18'),
                OFile(name='Sample31.gene_level.txt', size=507390, hash='293752c996a08ccac35146adaac9894540ae5ec1'),
                OFile(name='Sample31_hisens.cncf.txt', size=5124, hash='dce72617f478ebcb1e4644ed7e60e5282c303c1e'),
                OFile(name='Sample31_hisens.rds', size=625488, hash='5b88eb5926dd6968d56998c2160a104dad1bba65'),
                OFile(name='Sample31_hisens.seg', size=1920, hash='2845380bf9867e48d533dadf677f90d1cdc9b59c'),
                OFile(name='Sample31_purity.rds', size=625054, hash='2dfefbd3d975a7401bb1c3d7cc41de229e3c39fb'),
                OFile(name='Sample31_purity.seg', size=1507, hash='2a41207f14870350fed5cf50ad8ab6bf73e18c0f'),
                OFile(name='Sample31.qc.txt', size=1294, hash='a14a3d3a9a5e21a01af6969073f8996ed1147a2f')]),
            ODir(name='Sample14.s_C_Patient14_N001_d', items=[
                OFile(name='Sample14.s_C_Patient14_N001_d_hisens.ccf.portal.maf', size=49410903, hash='efe82e46d9d479491928630f3b1b5e656558ddb8'),
                OFile(name='Sample14.arm_level.txt', size=1641, hash='d850d38408368e609616d5201fa5475a385f40e2'),
                OFile(name='Sample14.txt', size=637, hash='b91d8eaa803dd642451c93008a678555d9885ca3'),
                OFile(name='Sample14.gene_level.txt', size=625826, hash='d94f2ded80d10e7c893f1b1b1153be5cc6c0a812'),
                OFile(name='Sample14_hisens.cncf.txt', size=9943, hash='4eeda7b2441d288083facd3c59c116f37d856712'),
                OFile(name='Sample14_hisens.rds', size=699789, hash='9ae2aca001adad670d21dd888966299393adce88'),
                OFile(name='Sample14_hisens.seg', size=3302, hash='103adc4518995aa4ffe3adcd198c51bb1ca890d1'),
                OFile(name='Sample14_purity.rds', size=699072, hash='adce7012ccade84cfb184f77afdfb433e6ee4e52'),
                OFile(name='Sample14_purity.seg', size=2529, hash='7a1d26f06fe1c9a5a39d503aafe0bb1949199194'),
                OFile(name='Sample14.qc.txt', size=1380, hash='812844a349c1b367e5130cde5c7f0a9243a6958c')]),
            ODir(name='Sample26.s_C_Patient26_N001_d', items=[
                OFile(name='Sample26.s_C_Patient26_N001_d_hisens.ccf.portal.maf', size=120843666, hash='bf397294a3be13cb7ff446f8c2e72afdda3d35bc'),
                OFile(name='Sample26.arm_level.txt', size=1906, hash='fa62a1159122112a40497458bdeb8ac627df6e55'),
                OFile(name='Sample26.txt', size=585, hash='86e3657fcb176aa921112acddedcf2de2a2e85aa'),
                OFile(name='Sample26.gene_level.txt', size=476434, hash='2ed5a12f1865f8a4e62733c979f420f11543a516'),
                OFile(name='Sample26_hisens.cncf.txt', size=4379, hash='592cf011d59919ca4fba6cb16fdb269ef704d547'),
                OFile(name='Sample26_hisens.rds', size=561304, hash='27dbf164fd45eaa9e4968c150e6f44ff87106e3a'),
                OFile(name='Sample26_hisens.seg', size=1672, hash='448494b4a2689012e678efd64cec78d2c65bf4b3'),
                OFile(name='Sample26_purity.rds', size=560910, hash='0dd7910e84089aa19e01842487a157e17f2fe4c1'),
                OFile(name='Sample26_purity.seg', size=1252, hash='a6e9097251da6efeeb6a4a31a4f3f0eb6c979ae2'),
                OFile(name='Sample26.qc.txt', size=1258, hash='acc344fda2bcb92ab5c0614963dde3b7abd1fc87')]),
            ODir(name='Sample8.s_C_Patient8_N001_d', items=[
                OFile(name='Sample8.s_C_Patient8_N001_d_hisens.ccf.portal.maf', size=52973498, hash='676990836edd84654beeb787976ce389bb1390ca'),
                OFile(name='Sample8.arm_level.txt', size=1859, hash='d152cdb6ee0639d97f6edda5d65025c96672d421'),
                OFile(name='Sample8.txt', size=575, hash='26b9c9f4bac527a2323a0278f866fd62165d44ee'),
                OFile(name='Sample8.gene_level.txt', size=470374, hash='b587d36fa9062c5f7f8ad9dc5bc3ec25c12e6b7c'),
                OFile(name='Sample8_hisens.cncf.txt', size=4252, hash='027a45325cefe43407cb8bfdab8f0aa182117893'),
                OFile(name='Sample8_hisens.rds', size=601597, hash='0edbe24060c88f6a72c0ad36dec57edd16b91d0a'),
                OFile(name='Sample8_hisens.seg', size=1586, hash='974f70c4db7b75a473683ed8056cdc175996d902'),
                OFile(name='Sample8_purity.rds', size=601171, hash='abb4c4f2d72b6fc831a6581bab53fbdff52dbf4a'),
                OFile(name='Sample8_purity.seg', size=1277, hash='ab2264b66fb5447ee26ab543d05335d6015855b7'),
                OFile(name='Sample8.qc.txt', size=1258, hash='07ca844c89078e77891333315736421c2aa41a53')]),
            ODir(name='Sample9.s_C_Patient9_N001_d', items=[
                OFile(name='Sample9.s_C_Patient9_N001_d_hisens.ccf.portal.maf', size=63376417, hash='208ab256b029bf4e5620768580ba0243259fcc1b'),
                OFile(name='Sample9.arm_level.txt', size=1865, hash='df0d168ef71e9ced38fd796e9319efd808fb2ebe'),
                OFile(name='Sample9.txt', size=586, hash='52fa384c83f4fb8432aecff39dc394d36f0e4dde'),
                OFile(name='Sample9.gene_level.txt', size=548958, hash='e12b127fdacef762bba59f8b3e7788d84761f566'),
                OFile(name='Sample9_hisens.cncf.txt', size=5174, hash='b44afe8f3010cd8dadf0c6d788ba8f46f3389aef'),
                OFile(name='Sample9_hisens.rds', size=709484, hash='c36f4254fb522fa5ea577d8393d181fac7da2b35'),
                OFile(name='Sample9_hisens.seg', size=1963, hash='60abece493d552bc0546b742192e257e84cc2b6c'),
                OFile(name='Sample9_purity.rds', size=708785, hash='9082db29d6508ff57fc7f90df57f9a660ed77a57'),
                OFile(name='Sample9_purity.seg', size=1346, hash='26ea1345a006b1aa61db8879639945b189c58dc4'),
                OFile(name='Sample9.qc.txt', size=1276, hash='27ed22c0f8031d15c3f01f4dc7dfb24da956b50d')]),
            ODir(name='Sample23.s_C_Patient23_N001_d', items=[
                OFile(name='Sample23.s_C_Patient23_N001_d_hisens.ccf.portal.maf', size=58240550, hash='9731fa1deb05048d75b1542960d3bf02d309a77d'),
                OFile(name='Sample23.arm_level.txt', size=1808, hash='dc7fe7703924005a80802016a9050ccf68fa8789'),
                OFile(name='Sample23.txt', size=598, hash='f83d988af6aa6122e95d11f79ca0453663fc17a0'),
                OFile(name='Sample23.gene_level.txt', size=731082, hash='e32df14c4b6fd265e44a3de48c3a92709b3d6b02'),
                OFile(name='Sample23_hisens.cncf.txt', size=7513, hash='2b4d1c69c7fd3ea87e1db120cf7bdfb5c4c63dbd'),
                OFile(name='Sample23_hisens.rds', size=690465, hash='d86cc62afce90765de03c326959aa2f7d0471dbd'),
                OFile(name='Sample23_hisens.seg', size=2701, hash='174b5943fc39481e511914a02b6ea61fbfb846a1'),
                OFile(name='Sample23_purity.rds', size=689404, hash='6e76701a8a130626c51fcb735eb7bb52899efaa9'),
                OFile(name='Sample23_purity.seg', size=1564, hash='1149af2192e3533a680f526ff029430392dbf71f'),
                OFile(name='Sample23.qc.txt', size=1361, hash='fcdc475e00780e231fe1974330d68ce058fcf27a')]),
            ODir(name='Sample6.s_C_Patient6_N001_d', items=[
                OFile(name='Sample6.s_C_Patient6_N001_d_hisens.ccf.portal.maf', size=86941095, hash='a0af34370cbd4845f7a14186d101cc7513d7e902'),
                OFile(name='Sample6.arm_level.txt', size=1482, hash='c9c9fceaca695d4e1352a17687296332468ebdf1'),
                OFile(name='Sample6.txt', size=794, hash='2c21f4e275ae5a8ce7117a1db19251b1a51f3419'),
                OFile(name='Sample6.gene_level.txt', size=681161, hash='5074e21603d5f430dbe6d620d61e849dfceeb8c5'),
                OFile(name='Sample6_hisens.cncf.txt', size=8656, hash='e409eac5197af981f9e809b22b87438cafff24f8'),
                OFile(name='Sample6_hisens.rds', size=681939, hash='5b62ebc56b9f7bcbdb6233cce2644a019bf48699'),
                OFile(name='Sample6_hisens.seg', size=2846, hash='4b44d3d52005f17c1df94486da271ce1f89f55a9'),
                OFile(name='Sample6_purity.rds', size=681844, hash='3b8fe0249f4d88df303e440249ae80286d256d60'),
                OFile(name='Sample6_purity.seg', size=2592, hash='b5dc2fa99d5a413c66e4d43a92ef2b4eb2e3480b'),
                OFile(name='Sample6.qc.txt', size=1358, hash='49629d1ebd5fcfa0cef9316e4420d2192a5bc829')]),
            ODir(name='Sample29.s_C_Patient29_N001_d', items=[
                OFile(name='Sample29.s_C_Patient29_N001_d_hisens.ccf.portal.maf', size=86674879, hash='d34d89c2c15bfec282b5465235162c2304d40cab'),
                OFile(name='Sample29.arm_level.txt', size=1680, hash='231f860c7a48cdb039a6a427ac6bf533550903d9'),
                OFile(name='Sample29.txt', size=645, hash='cd658ef98ff6a036afcaa6ec0d90b6bb12562bd0'),
                OFile(name='Sample29.gene_level.txt', size=663879, hash='f32cb0da212dc0f909e637b279d86a858580447d'),
                OFile(name='Sample29_hisens.cncf.txt', size=11344, hash='9c7f080b7b2424e97815e0b0450ecee7d9290763'),
                OFile(name='Sample29_hisens.rds', size=648985, hash='b7612862fad77f05bc29dc68b9146e972677589c'),
                OFile(name='Sample29_hisens.seg', size=3995, hash='df47d69ae42e1ad2a7b3c265e00f0d4a321a9724'),
                OFile(name='Sample29_purity.rds', size=646477, hash='b7fe6d04c1a6312f69152a474a2c5726742f11e4'),
                OFile(name='Sample29_purity.seg', size=1757, hash='8452f272cb686344bf006457731396cdfeb26fa0'),
                OFile(name='Sample29.qc.txt', size=1373, hash='ff014829c3cd5f4de5b1139320b92fe0f5bf8ce1')]),
            ODir(name='Sample4.s_C_Patient4_N001_d', items=[
                OFile(name='Sample4.s_C_Patient4_N001_d_hisens.ccf.portal.maf', size=49423425, hash='0ac006f61cc19509e837ed780eba5ad07621d1d0'),
                OFile(name='Sample4.arm_level.txt', size=1922, hash='2223703c91c4233f68b9b645eb8d169a8b495bbf'),
                OFile(name='Sample4.txt', size=614, hash='31d946ef4e803893b469d74271c1d9f3445699c7'),
                OFile(name='Sample4.gene_level.txt', size=809724, hash='2a81fe2b40c5fadcada16003d12dca2f7296ac29'),
                OFile(name='Sample4_hisens.cncf.txt', size=8255, hash='2869f7abe4721e893c1e89a9dc5a995cbdaff259'),
                OFile(name='Sample4_hisens.rds', size=837983, hash='9403a28138d9571f0f3d1957d0db6df8fbbb19f4'),
                OFile(name='Sample4_hisens.seg', size=2863, hash='c787e0dbc9e3692b815a8af63e3ee22be1f041ce'),
                OFile(name='Sample4_purity.rds', size=836693, hash='3e4bba214a9ad5d0a1aa8a164382adb78a73bdfd'),
                OFile(name='Sample4_purity.seg', size=1681, hash='c96476425746d373de5b44322f1708e6ac9c1129'),
                OFile(name='Sample4.qc.txt', size=1336, hash='47e06f5834afafa614851e0e55e1b12fabe12914')]),
            ODir(name='Sample30.s_C_Patient30_N001_d', items=[
                OFile(name='Sample30.s_C_Patient30_N001_d_hisens.ccf.portal.maf', size=94320907, hash='dab1ee980963e8cf069909488c03bd0c6d15c9aa'),
                OFile(name='Sample30.arm_level.txt', size=1894, hash='ba5b095b64e8c7ff422719ad8807c9a49db023ee'),
                OFile(name='Sample30.txt', size=601, hash='9d0e82d5919e851c9470ff434cebf8d1ab1fefd0'),
                OFile(name='Sample30.gene_level.txt', size=790775, hash='e500d61c1f2a5b531788f8e40104206a998f5908'),
                OFile(name='Sample30_hisens.cncf.txt', size=8523, hash='ea38582b14a8924b98e64a8c2b0237afd8552b63'),
                OFile(name='Sample30_hisens.rds', size=768241, hash='77b8e39b30028cd6a9681793edfe24755159f724'),
                OFile(name='Sample30_hisens.seg', size=3283, hash='8d1c7400139769eaed3620485a582fe501dec98c'),
                OFile(name='Sample30_purity.rds', size=765993, hash='a6aac14b363765f613bfb72eb642ed96f1a6ba5f'),
                OFile(name='Sample30_purity.seg', size=1362, hash='49c0f716e12a043e0c31d763b6b5de49960b9ef6'),
                OFile(name='Sample30.qc.txt', size=1344, hash='484df5e7b9561a3f477043078bffe84f38b0eb55')]),
            ODir(name='Sample27.s_C_Patient27_N001_d', items=[
                OFile(name='Sample27.s_C_Patient27_N001_d_hisens.ccf.portal.maf', size=45005876, hash='8784e9e95c6a392d92ba9bbf71ed87462132496f'),
                OFile(name='Sample27.arm_level.txt', size=1906, hash='02d2da95e89bb017e1603c484ea7b08b9b9fb856'),
                OFile(name='Sample27.txt', size=582, hash='83397779b3d37b2177b0e573af9bcb2c3d02c6d5'),
                OFile(name='Sample27.gene_level.txt', size=445227, hash='4792260db0f74a843aadd7578e1cc8753082ed3f'),
                OFile(name='Sample27_hisens.cncf.txt', size=3406, hash='4a8df4aece8adbeaa043c060ed5ad952a0c95ef1'),
                OFile(name='Sample27_hisens.rds', size=599275, hash='270aa7723ccf97cbc68e02010713b0dea8ccff6f'),
                OFile(name='Sample27_hisens.seg', size=1257, hash='a9eef526748e2f1cb8a6cc634761c2b5bddca704'),
                OFile(name='Sample27_purity.rds', size=599333, hash='ffb261d90559d1d656f5bd0711f487d6f427257c'),
                OFile(name='Sample27_purity.seg', size=1257, hash='a9eef526748e2f1cb8a6cc634761c2b5bddca704'),
                OFile(name='Sample27.qc.txt', size=1268, hash='d3f5a4522297d4dfa7503ecb9764aae05d8df418')]),
            ODir(name='Sample2.s_C_Patient2_N001_d', items=[
                OFile(name='Sample2.s_C_Patient2_N001_d_hisens.ccf.portal.maf', size=78410554, hash='854608265cb3f020a3000a54074cb15727dbe885'),
                OFile(name='Sample2.arm_level.txt', size=1742, hash='050255a51a1476f211a71acbbb52f73c2541dae2'),
                OFile(name='Sample2.txt', size=623, hash='128de962ad52ccbf06f975b3039f0b0c54d417d2'),
                OFile(name='Sample2.gene_level.txt', size=655694, hash='45318af0025bfbd1a8cb97628f9729f44e0ccc85'),
                OFile(name='Sample2_hisens.cncf.txt', size=7658, hash='a854b06c09f30966bc9137d6a14d0860c6138ccf'),
                OFile(name='Sample2_hisens.rds', size=682739, hash='327ca88e5bfb7f135de4dfc60e06cf25b7967c26'),
                OFile(name='Sample2_hisens.seg', size=2644, hash='6275921e635ec4044544d6c9343c67aa7dcd142c'),
                OFile(name='Sample2_purity.rds', size=682025, hash='6261f146e6cad666f38862ed99ab0bbcbe1aab2a'),
                OFile(name='Sample2_purity.seg', size=1928, hash='2e62b02dd63e4f1e6193480f3f741343ded1cb37'),
                OFile(name='Sample2.qc.txt', size=1335, hash='fd593d3a1babf3903c95875114fba069812b907c')]),
            ODir(name='Sample22.s_C_Patient22_N001_d', items=[
                OFile(name='Sample22.s_C_Patient22_N001_d_hisens.ccf.portal.maf', size=232656074, hash='f9230b388b7688c257247b55a6d07e2a3cbe357b'),
                OFile(name='Sample22.arm_level.txt', size=1877, hash='0d7680007e0c1d58efc39b6f3d1c4e907a09243e'),
                OFile(name='Sample22.txt', size=600, hash='b43007c41f7810d5c0ab8cd3ccb05c204f3e6a76'),
                OFile(name='Sample22.gene_level.txt', size=578835, hash='0cfbcd2c032aaf9ae8288c6c9c08b2dc9b9852d7'),
                OFile(name='Sample22_hisens.cncf.txt', size=3398, hash='c4c838462fd1c806cab5d44a60a7f5f1f97c1e9d'),
                OFile(name='Sample22_hisens.rds', size=616983, hash='8fadbe6da539a385d810036c7bf3f3d50b7cc071'),
                OFile(name='Sample22_hisens.seg', size=1256, hash='673c77974d839ccc3464946cb52b5e8726b06a6f'),
                OFile(name='Sample22_purity.rds', size=617039, hash='ebe499a84bd8f7425746ed4c17f29260b6a18718'),
                OFile(name='Sample22_purity.seg', size=1256, hash='673c77974d839ccc3464946cb52b5e8726b06a6f'),
                OFile(name='Sample22.qc.txt', size=1334, hash='9dffb6ce9b710bcd25f664e47af54b532644b00e')]),
            ODir(name='Sample25.s_C_Patient25_N001_d', items=[
                OFile(name='Sample25.s_C_Patient25_N001_d_hisens.ccf.portal.maf', size=89648628, hash='6006cc0fdc31821292d9dca5281128844ca4b58b'),
                OFile(name='Sample25.arm_level.txt', size=1852, hash='9cb1c4584de1dfe10411a481bdea867db00d8582'),
                OFile(name='Sample25.txt', size=597, hash='36c9957942bd99e25098c533cc09102c76000c60'),
                OFile(name='Sample25.gene_level.txt', size=481755, hash='809c7796a101c0fb4e993a9b3c7d4f84ac8e8a15'),
                OFile(name='Sample25_hisens.cncf.txt', size=3852, hash='120dbd0b92153c30abdd46ae5ef353c8c1306a77'),
                OFile(name='Sample25_hisens.rds', size=583358, hash='b56336324db262f72cdbb5968ef272da68392df3'),
                OFile(name='Sample25_hisens.seg', size=1403, hash='184f5c5f23600a2ed7d999609cac8b333a9a0587'),
                OFile(name='Sample25_purity.rds', size=583275, hash='7bb7b780313cf1fe7126f3246c2f51c998e224f0'),
                OFile(name='Sample25_purity.seg', size=1302, hash='6b6d5c266e47c84cab638f4c8345fe9c78ec5825'),
                OFile(name='Sample25.qc.txt', size=1328, hash='0613673c776fc4ede4a63e33c9e5193e1c603fa7')])],
            dir=output_dir),
        'facets_failed_pairs': [],
        'portal_dir': ODir(name='portal', items=[
            OFile(name='meta_clinical_sample.txt', size=135, hash='eab82493f207bd3b7ae6a71bec253ede9371e160'),
            OFile(name='data_clinical_patient.txt', size=492, hash='056fe0afc68169191c0981faa9c103a73d0cd2c8'),
            OFile(name='data_clinical_sample.txt', size=6397, hash='99309416dfb2aa14818e9dbda138d2c4a035c829'),
            OFile(name='meta_study.txt', size=138, hash='3103ac7f24c557ec35ef3692e3cdb166427a1590'),
            OFile(name='meta_clinical_patient.txt', size=137, hash='7c7ba5be461f46f9b4728025717beb4b4ceb9482'),
            OFile(name='meta_CNA.txt', size=265, hash='53963142b01fb5e1ae2c92feeff39e0027559fb1'),
            OFile(name='meta_fusions.txt', size=222, hash='fec37f17f5adc2a5861b3c382f5f01d810c67145'),
            OFile(name='meta_mutations_extended.txt', size=265, hash='cbd9bf1d92c587887f9f823ad799de00ce28b152'),
            OFile(name='11089_G_meta_cna_hg19_seg.txt', size=190, hash='03acd909d48042090e5e941335b4756b30133841'),
            OFile(name='data_CNA.txt', size=140732, hash='5b0f504d0fe9117841407864e3167008f5713e7e'),
            OFile(name='data_CNA.ascna.txt', size=227182, hash='61b435a25a2df0e47b65709e197757da1662f67d'),
            OFile(name='data_mutations_extended.txt', size=49008, hash='6cffa32a45655147c86de5aacd95a4f9f807cbc0'),
            OFile(name='11089_G_data_cna_hg19.seg', size=69171, hash='7d13e1422b9d7112e60f9ab748d524d7446d1dde'),
            OFile(name='data_fusions.txt', size=406, hash='4c9e29ead2ddd262dc3877dcc563c6c1d5c624ad'),
            ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=448, hash='da1d6dce1a3d9694bfb74bb99b1cda1766531164'),
                OFile(name='cases_cnaseq.txt', size=528, hash='0bf89afc7d8409185288c8bdee45ffb07daa4b08'),
                OFile(name='cases_cna.txt', size=460, hash='dfc9325922ceef49a427a61362366507b6331948'),
                OFile(name='cases_sequenced.txt', size=473, hash='6aa2ba1b519020ada95acb77b0900be172f306ff')]),
                OFile(name='report.html', size=1026290, hash='8c70ab473d323df777c86d330d6e064eaf580cd6')],
            dir=output_dir)
        }
        self.maxDiff = None
        self.assertCWLDictEqual(output_json, expected_output)
        self.assertNumMutations(os.path.join(output_dir, 'analysis', '11089_G.muts.maf'), 34)
        self.assertNumMutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'), 27)
        self.assertSampleValues(
            os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
            value_fieldname = "CMO_TMB_SCORE",
            expected_values = {
            'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': '5.5', 'Sample1': '47.5'
            })



