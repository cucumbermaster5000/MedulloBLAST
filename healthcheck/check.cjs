const {chromium}=require('playwright');
const fs=require('fs');const path=require('path');
const out=process.env.HEALTH_OUTPUT||'health-results';fs.mkdirSync(out,{recursive:true});
const gene='HLX',url='https://cucumbermaster5000-medulloblast.share.connect.posit.cloud/';
(async()=>{let browser,page;const result={started:new Date().toISOString(),url,gene,checks:[],errors:[]};const record=(name,detail)=>{result.checks.push({name,detail});console.log('PASS '+name+': '+detail)};
try{
 browser=await chromium.launch(process.env.HEALTH_BROWSER_CHANNEL?{headless:true,channel:process.env.HEALTH_BROWSER_CHANNEL}:{headless:true});page=await browser.newPage({acceptDownloads:true,viewport:{width:1440,height:960}});page.setDefaultTimeout(180000);page.on('pageerror',e=>result.errors.push(e.message));
 const response=await page.goto(url,{waitUntil:'domcontentloaded',timeout:120000});if(!response.ok())throw Error('App HTTP '+response.status());
 await page.waitForFunction(()=>window.Shiny?.shinyapp?.$socket?.readyState===1,null,{timeout:120000});record('connect','HTTP '+response.status());
 await page.locator('#gene').fill(gene);await page.locator('#run').click();
 await page.locator('#status .analysis-ready').waitFor({timeout:420000});let status=await page.locator('#status .analysis-ready').innerText();if(status.includes('Unavailable:')||!status.includes(gene))throw Error('Analysis status: '+status);record('loading',status);
 async function tab(name){await page.getByRole('tab',{name,exact:true}).click()}
 async function text(id,pattern){await page.waitForFunction(({id,pattern})=>new RegExp(pattern,'i').test(document.getElementById(id)?.textContent||''),{id,pattern},{timeout:180000});return (await page.locator('#'+id).innerText()).slice(0,250)}
 async function plot(id){await page.waitForFunction(id=>{let i=document.querySelector('#'+id+' img');return i&&i.complete&&i.naturalWidth>0},id,{timeout:180000})}
 await tab('MB Subgroups');await text('plot_heading',gene+' expression');await plot('subgroup_plot');record('Cavalli','subgroup plot');
 await tab('MB Subtypes');await plot('subtype_plot');record('subtypes','plot rendered');
 await tab('Pediatric Pan-Cancer');record('Pfister',await text('pfister_summary',gene+'.*Pfister'));await plot('pfister_plot');
 await tab('Survival');await page.waitForFunction(()=>{const t=document.body.innerText;return t.includes('Verified usable:')||t.includes('could not be verified')},null,{timeout:180000});const survival=await page.locator('.tab-pane.active:visible').first().innerText();if(!survival.includes('Verified usable:')||survival.includes('could not be verified'))throw Error('Survival metadata not verified: '+survival.slice(0,800));await page.waitForFunction(()=>Array.from(document.querySelectorAll('.tab-pane.active img')).some(i=>i.complete&&i.naturalWidth>0),null,{timeout:180000});record('survival','verified metadata and plot');
 await tab('snRNA-seq Developing Cerebellum');record('cerebellum human',await text('cerebellum_expression_state_human','expressing /'));record('cerebellum mouse',await text('cerebellum_expression_state_mouse','expressing /'));
 await tab('DepMap');record('DepMap',await text('depmap_status','DepMap loaded for '+gene));await plot('depmap_dependency_plot');
 await tab('Publications');record('PubMed',await text('publication_status',gene+' - .*matching PubMed results'));
 await tab('Functional Biology');record('functional biology',await text('biology_status','loaded for '+gene));
 await tab('ProteinAtlas.org Profile');record('Protein Atlas',await text('hpa_status',gene+' - Human Protein Atlas'));
 await tab('MB Subgroups');const pending=page.waitForEvent('download',{timeout:180000});await page.locator('#cavalli_subgroups_csv').click();const download=await pending;const file=path.join(out,'HLX_Cavalli.csv');await download.saveAs(file);const csv=fs.readFileSync(file,'utf8');const lines=csv.trim().split(/\r?\n/);if(!download.suggestedFilename().startsWith('HLX_')||lines.length!==764||!lines[0].includes('gene')||!lines[1].includes('HLX'))throw Error('Unexpected Cavalli CSV content');record('CSV download',`${lines.length-1} patient rows`);
 if(result.errors.length)throw Error('Browser errors: '+result.errors.join('; '));result.passed=true;
}catch(e){result.passed=false;result.errors.push(String(e));if(page){try{await page.screenshot({path:path.join(out,'failure.png'),fullPage:true,timeout:15000});fs.writeFileSync(path.join(out,'page.txt'),await page.locator('body').innerText({timeout:15000}));}catch{}}console.error(e);process.exitCode=1;
}finally{result.finished=new Date().toISOString();fs.writeFileSync(path.join(out,'result.json'),JSON.stringify(result,null,2));if(browser)await browser.close()}
})();
