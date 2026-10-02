// Only real user interaction resets the server-side inactivity deadline.
(function () {
  let connected=false, expired=false, serial=0, lastSent=0, pending;
  function sendActivity() {
    pending=undefined;
    if(!connected||expired||!window.Shiny)return;
    lastSent=Date.now();
    Shiny.setInputValue('mb_user_activity',++serial,{priority:'event'});
  }
  function activity(event) {
    if(!event.isTrusted||document.visibilityState==='hidden'||expired)return;
    if(Date.now()-lastSent>=1000)sendActivity();
    else if(!pending)pending=setTimeout(sendActivity,1000-(Date.now()-lastSent));
  }
  ['pointerdown','pointermove','keydown','wheel','touchstart','scroll'].forEach(name=>
    document.addEventListener(name,activity,{capture:true,passive:true}));
  $(document).on('shiny:connected',function(){connected=true;sendActivity();});
  $(document).on('shiny:disconnected',function(){connected=false;clearTimeout(pending);});
  Shiny.addCustomMessageHandler('mb-session-policy',function(message){
    const box=document.getElementById('mb-idle-warning');
    const text=document.getElementById('mb-idle-warning-text');
    if(message.type==='active'){if(box)box.style.display='none';return;}
    if(message.type==='warning'){
      if(text)text.textContent='Your inactive session will disconnect in about '+message.seconds+' seconds.';
      if(box)box.style.display='block';return;
    }
    if(message.type==='expired'){
      expired=true;clearTimeout(pending);
      const panel=document.createElement('div');panel.id='mb-session-expired';panel.setAttribute('role','alertdialog');
      panel.setAttribute('aria-label','Session disconnected');
      panel.style.cssText='position:fixed;inset:0;z-index:2147483647;background:rgba(255,255,255,.97);display:grid;place-content:center;padding:32px;text-align:center;color:#29384B';
      const title=document.createElement('h2');title.textContent='Session disconnected due to inactivity';
      const detail=document.createElement('p');detail.textContent='Your session was closed to free capacity for other users. Reconnect to run a new analysis.';
      const button=document.createElement('button');button.className='btn btn-primary';button.textContent='Reconnect';button.addEventListener('click',()=>location.reload());
      panel.append(title,detail,button);document.body.append(panel);
    }
  });
})();
